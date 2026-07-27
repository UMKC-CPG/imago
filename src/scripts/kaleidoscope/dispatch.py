## SPDX-License-Identifier: ECL-2.0
## Copyright (c) 2026 Paul Rulis

"""kaleidoscope.dispatch -- the flight driver
(DESIGN 6.2.3; PSEUDOCODE 13.5).

``dispatch`` walks the units, consults the cache, dispatches
the misses through an executor, and gathers the results with
per-future exception capture so a single failure never aborts
the flight (VISION Principle 10).  Resuming a flight is just
re-running it: the cache hit-test skips the units already
``done`` and re-dispatches the rest.

Two executors realize the dispatch:

- ``LocalExecutor`` runs units synchronously, in process -- no
  parallelism, no Parsl.  It is the default when a flight
  carries no ``parsl_config``, and is what tests and Parsl-less
  environments use.
- ``ParslExecutor`` dispatches units as Parsl ``python_app``\\s
  (DESIGN 6.2.3) and is selected when the flight supplies a
  Parsl ``Config``.  Choosing the executor by the presence of
  ``parsl_config`` keeps Parsl central (VISION Goal 4) without
  making it a hard import-time dependency.

Dispatch is two phases, and both are public (DESIGN 6.2.3): a
``send_off`` that launches a chosen set of units and returns one
future per unit without waiting, and a ``collect`` that resolves a
single future into its terminal status and report entry.  ``dispatch``
itself is just ``send_off`` followed by collecting every future in
unit order -- the convenience form for a one-shot fan-out.  A
control-loop client (the k-point climb) instead sends the rungs it
has decided and uses ``collect_next`` to take whichever lands first
and send that chain's successor at once, never waiting out a whole
batch.  ``collect_next`` is domain-ignorant -- it polls each future's
``done()`` and knows only futures, never k-points -- so the choice of
what to send next stays in the client (Principle 12).
"""

import os
import time
from datetime import datetime

from .model import FlightReport, ReportEntry
from .workspace import (unit_run_dir, validate_flight,
                        serialize_flight, write_status,
                        read_status)
from .cache import is_cache_hit, write_cache_key
from .wingbeats import resolve_wingbeat


# How long ``collect_next`` sleeps between scans of the outstanding
#   futures when none has landed yet (seconds).  The local executor
#   finishes synchronously, so its futures are already done and this
#   sleep is never reached; only a real cluster wait pays it, where a
#   half-second against multi-second-to-minutes rungs is negligible.
_COLLECT_POLL_SECONDS = 0.5


class TaskLost(Exception):
    """A unit whose executor task vanished with no WingbeatOutcome --
    a cluster-side loss (manager/worker death, expired
    allocation).  Mapped to the ``lost`` status, which is
    distinct from ``failed`` (a unit that ran and reported, or
    whose wingbeat raised) (DESIGN 6.2.4)."""
    pass


def now_iso():
    """Second-resolution ISO timestamp for status.toml fields."""
    return datetime.now().isoformat(timespec="seconds")


def _is_lost(exc):
    """Heuristic, version-robust classification of an executor
    exception as a cluster-side loss.  Parsl signals a vanished
    task with exception types whose names mention loss (e.g.
    ``ManagerLost``, ``WorkerLost``) or a block in a bad state
    (``BadStateException``).  Everything else is an ordinary
    failure."""
    name = type(exc).__name__
    return ("Lost" in name) or (name == "BadStateException")


# ------------------------------------------------------------------
#  The unit task (module-level so a worker process can import it)
# ------------------------------------------------------------------

def _execute_wingbeat_task(unit, wingbeat_dir, default_wingbeat):
    """Execute one unit and return its ``WingbeatOutcome``.  Runs on a
    worker (a thread, or a remote process under a multi-process
    executor), so it writes the ``running`` and the terminal
    status into the run directory itself.  Resolving the wingbeat
    by name here is what lets a worker process reconstruct it
    after import (DESIGN 6.2.2, 6.2.3)."""
    wingbeat_name = unit.wingbeat or default_wingbeat
    write_status(wingbeat_dir, id=unit.id, calc=unit.calc,
                 status="running", wingbeat=wingbeat_name,
                 started_at=now_iso())
    wingbeat = resolve_wingbeat(wingbeat_name)
    outcome = wingbeat.run(unit, wingbeat_dir)
    write_status(wingbeat_dir, id=unit.id, calc=unit.calc,
                 status=("done" if outcome.ok else "failed"),
                 detail=outcome.detail, finished_at=now_iso(),
                 runtime_seconds=outcome.runtime_seconds,
                 message=outcome.message)
    return outcome


# ------------------------------------------------------------------
#  Executors
# ------------------------------------------------------------------

class _LocalFuture:
    """A trivial future for the synchronous executor: it holds a
    value or an error and re-raises the error on ``result()``,
    mirroring the contract of a real future -- ``result()`` and
    ``done()``.  With no value and no error it is the already-done
    placeholder a cache hit returns (:func:`completed_future`)."""

    def __init__(self, value=None, error=None):
        self._value = value
        self._error = error

    def result(self):
        if self._error is not None:
            raise self._error
        return self._value

    def done(self):
        # The work already ran synchronously (or this is a cache-hit
        #   placeholder), so a local future is born finished.  This
        #   is what lets a local run never reach collect_next's poll
        #   sleep -- its futures are always ready on the first scan.
        return True


class LocalExecutor:
    """Run units synchronously in the current process (DESIGN
    6.2.3).  ``submit_unit`` runs the task immediately and
    captures any exception so it surfaces on ``result()`` exactly
    as a real future would -- letting the dispatcher's
    complete-and-report logic stay identical across executors."""

    def submit_unit(self, unit, wingbeat_dir, default_wingbeat):
        try:
            value = _execute_wingbeat_task(unit, wingbeat_dir, default_wingbeat)
            return _LocalFuture(value=value)
        except Exception as err:          # noqa: BLE001
            return _LocalFuture(error=err)

    def close(self):
        pass


class _ParslFuture:
    """Adapter over a Parsl ``AppFuture`` that converts a
    cluster-side loss into ``TaskLost`` while passing every other
    exception through unchanged."""

    def __init__(self, app_future):
        self._app_future = app_future

    def result(self):
        try:
            return self._app_future.result()
        except Exception as err:          # noqa: BLE001
            if _is_lost(err):
                raise TaskLost(str(err)) from err
            raise

    def done(self):
        # Delegate to the underlying Parsl AppFuture (a
        #   concurrent.futures.Future), so a caller can poll for
        #   completion without blocking on result() -- this is what
        #   collect_next needs to find whichever rung landed first.
        return self._app_future.done()


class ParslExecutor:
    """Dispatch units as Parsl ``python_app``\\s (DESIGN 6.2.3).
    Parsl is imported lazily so the package works without it;
    constructing this executor requires Parsl to be installed and
    a ``Config`` to load.  The same code serves a laptop
    (a ThreadPoolExecutor config) and a cluster (a
    HighThroughputExecutor + SLURM provider config) -- only the
    ``Config`` changes."""

    def __init__(self, parsl_config):
        import parsl
        from parsl import python_app
        self._parsl = parsl
        parsl.load(parsl_config)
        # Wrap the module-level task as a Parsl app once; calling
        #   the app returns an AppFuture per unit.
        self._app = python_app(_execute_wingbeat_task)

    def submit_unit(self, unit, wingbeat_dir, default_wingbeat):
        return _ParslFuture(
            self._app(unit, wingbeat_dir, default_wingbeat)
        )

    def close(self):
        # Clean up the data-flow kernel so a later flight can
        #   load a fresh config in this process.
        self._parsl.dfk().cleanup()
        self._parsl.clear()


def make_executor(parsl_config):
    """Turn a resolved dispatch config into a live executor
    (PSEUDOCODE 13.5).  This is the single place the choice lives,
    so the driver and any caller that pins its own executor stay in
    agreement about what a given config means.

    A cluster ``parsl_config`` yields a :class:`ParslExecutor`;
    constructing it loads Parsl, which starts one coordinator
    process and its pool of SLURM workers (DESIGN 6.2.11).  The
    local opt-out -- ``parsl_config`` is ``None`` -- yields an
    in-process :class:`LocalExecutor` that touches no scheduler.

    A client that dispatches MANY flights under one config must
    call this ONCE and pin the result to every :func:`dispatch`.
    The producer's climb is exactly that shape -- a pre-flight
    batch and then one flight per continuation round (DESIGN
    3.12.5) -- and a Parsl config's executor is single-use: once
    closed, its coordinator and worker pool are gone, so handing
    the same config to a second dispatch would try to restart a
    torn-down pool.  Building the executor here and sharing it lets
    the whole run ride one warm pool.  A one-shot flight instead
    lets :func:`dispatch` build and close its own.
    """
    if parsl_config is not None:
        return ParslExecutor(parsl_config)
    return LocalExecutor()


# ------------------------------------------------------------------
#  Dispatch helpers
# ------------------------------------------------------------------

def _prepare_miss(flight, unit, wingbeat_dir):
    """Set up a cache miss for launch: create the run directory,
    snapshot the cache key, and mark the unit ``queued`` (DESIGN
    6.2.5)."""
    os.makedirs(wingbeat_dir, exist_ok=True)
    write_cache_key(wingbeat_dir, unit)
    write_status(wingbeat_dir, id=unit.id, calc=unit.calc,
                 status="queued",
                 wingbeat=(unit.wingbeat or flight.default_wingbeat),
                 submitted_at=now_iso())


def report_entry_from_status(unit, wingbeat_dir):
    """Build a ReportEntry from the run directory's terminal
    ``status.toml`` -- the single source of truth (DESIGN
    6.2.6)."""
    status = read_status(wingbeat_dir) or {}
    return ReportEntry(
        id=unit.id, calc=unit.calc,
        status=status.get("status", "unknown"),
        detail=status.get("detail"),
        wingbeat_dir=wingbeat_dir,
        runtime_seconds=status.get("runtime_seconds"),
        message=status.get("message"),
    )


def completed_future():
    """An already-done future for a cache hit: no task was submitted,
    so ``result()`` is a no-op and ``done()`` is True.  Returning one
    lets a hit sit in the outstanding set exactly like a miss; the
    entry is rebuilt from the existing ``status.toml`` when the hit is
    collected (DESIGN 6.2.3)."""
    return _LocalFuture()


def dispatch_unit(flight, unit, executor, force):
    """Launch ONE unit and return its future (PSEUDOCODE 13.5).  A
    cache hit submits no task and returns an already-done future
    (:func:`completed_future`); a miss is prepared -- run directory,
    cache-key snapshot, ``queued`` status -- and handed to the
    executor.

    ``force`` bypasses the run-reuse cache so even a still-valid
    ``done`` unit re-launches (DESIGN 6.2.5).  It rides here, on the
    driver, because the cache it governs is the driver's -- not the
    executor's -- so it is independent of which executor runs the
    unit."""
    wingbeat_dir = unit_run_dir(flight, unit)
    if not force and is_cache_hit(unit, wingbeat_dir):
        return completed_future()
    _prepare_miss(flight, unit, wingbeat_dir)
    return executor.submit_unit(
        unit, wingbeat_dir, flight.default_wingbeat)


# ------------------------------------------------------------------
#  The flight driver -- two public phases, and the one-shot wrapper
# ------------------------------------------------------------------

def send_off(flight, units, executor, force=False):
    """Phase 1, callable on its own: launch ``units`` and return one
    ``(unit, future)`` pair per unit WITHOUT waiting on any of them
    (DESIGN 6.2.3).  ``units`` is the subset to launch now --
    ``flight.units`` for a one-shot fan-out, one climb round's newly
    decided rungs for the adaptive climb (PSEUDOCODE 4e.7).

    The flight's WHOLE unit list is (re)serialized here, not just the
    launched subset, so ``flight.toml`` records every unit asked for
    even as a climb's list grows across successive sends (DESIGN
    7.7)."""
    validate_flight(flight)
    os.makedirs(flight.root, exist_ok=True)
    serialize_flight(flight)
    outstanding = []                  # (unit, future)
    for unit in units:
        future = dispatch_unit(flight, unit, executor, force)
        outstanding.append((unit, future))
    return outstanding


def collect(flight, unit, future):
    """Phase 2 for a SINGLE unit: resolve its future, write its
    terminal status, build its report entry, and fire the optional
    per-unit outcome hook (DESIGN 6.2.3, 6.2.6).

    A ``TaskLost`` becomes the ``lost`` status; any other exception
    (including a wingbeat that raised on the worker) becomes
    ``failed``; a clean return leaves the terminal status the task
    itself already wrote (DESIGN 6.2.3, 6.2.4).  The outcome hook
    fires in LANDING order -- the order units are actually collected,
    not unit order -- so a control-loop consumer (the climb) sees each
    rung the moment it is collected."""
    wingbeat_dir = unit_run_dir(flight, unit)
    try:
        future.result()
    except TaskLost as lost:
        write_status(wingbeat_dir, id=unit.id, calc=unit.calc,
                     status="lost", finished_at=now_iso(),
                     message=str(lost) or "cluster-side loss")
    except Exception as err:              # noqa: BLE001
        write_status(wingbeat_dir, id=unit.id, calc=unit.calc,
                     status="failed", finished_at=now_iso(),
                     message=str(err))
    entry = report_entry_from_status(unit, wingbeat_dir)
    if flight.on_outcome is not None:
        flight.on_outcome(entry)
    return entry


def collect_next(flight, outstanding):
    """Wait for WHICHEVER outstanding rung lands first, collect it,
    and return ``(unit, entry, remaining)`` -- the landed unit, its
    report entry, and the outstanding list with that unit removed
    (DESIGN 6.2.3).

    Domain-ignorant: it polls each future's ``done()`` and knows only
    futures, never k-points, so a control-loop client (the climb)
    keeps the decision of what to send next to itself (Principle 12).
    A cache hit is an already-done future, so it is returned first
    with no wait.  The local executor's futures are always done, so
    the poll sleep below is reached only on a real cluster wait, where
    a fraction of a second against multi-second rungs is negligible."""
    while True:
        for index, (unit, future) in enumerate(outstanding):
            if future.done():
                entry = collect(flight, unit, future)
                remaining = (outstanding[:index]
                             + outstanding[index + 1:])
                return unit, entry, remaining
        time.sleep(_COLLECT_POLL_SECONDS)


def dispatch(flight, executor=None, force=False):
    """Run every unit in the flight and return a FlightReport
    (DESIGN 6.2.3): the one-shot convenience form of the two public
    phases -- send every unit off, then collect them all in unit
    order.  Behaviour is identical to the pre-split driver, so every
    existing caller is unchanged; a control-loop client (the climb)
    uses :func:`send_off` and :func:`collect_next` directly instead.

    When ``executor`` is None one is chosen from the flight: a
    ``ParslExecutor`` if it carries a ``parsl_config``, else a
    ``LocalExecutor``.  A caller may pass an executor explicitly
    (tests do, to pin the path; the climb does, to share one warm
    pool across the whole run).

    ``force`` bypasses the run-reuse cache (DESIGN 6.2.5): every unit
    re-launches even if a completed ``status.toml`` already exists.
    The switch lives here, on the driver, because the cache it governs
    is owned by the driver -- not by the executor and not by any one
    client."""
    owns_executor = executor is None
    if executor is None:
        executor = make_executor(flight.parsl_config)
    try:
        outstanding = send_off(flight, flight.units, executor, force)
        # Collect in send order, which is unit order; the outcome
        #   hook fires inside collect, once per entry.
        entries = [collect(flight, unit, future)
                   for unit, future in outstanding]
    finally:
        if owns_executor:
            executor.close()
    return FlightReport(entries=entries)
