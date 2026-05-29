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
"""

import os
from datetime import datetime

from .model import FlightReport, ReportEntry
from .workspace import (unit_run_dir, validate_flight,
                        serialize_flight, write_status,
                        read_status)
from .cache import is_cache_hit, write_cache_key
from .wingbeats import resolve_wingbeat


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
    mirroring the contract of a real future."""

    def __init__(self, value=None, error=None):
        self._value = value
        self._error = error

    def result(self):
        if self._error is not None:
            raise self._error
        return self._value


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


def _collect(flight, unit, future):
    """Resolve one future, recording its terminal status.  A
    ``TaskLost`` becomes the ``lost`` status; any other exception
    (including a wingbeat that raised on the worker) becomes
    ``failed``; a clean return leaves the terminal status the
    task itself already wrote (DESIGN 6.2.3, 6.2.4)."""
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
    return report_entry_from_status(unit, wingbeat_dir)


def _fire(flight, entry):
    """Invoke the optional per-unit outcome callback."""
    if flight.on_outcome is not None:
        flight.on_outcome(entry)


# ------------------------------------------------------------------
#  The flight driver
# ------------------------------------------------------------------

def dispatch(flight, executor=None):
    """Run every unit in the flight and return a FlightReport
    (DESIGN 6.2.3).  Cache hits are reported straight from their
    existing ``status.toml``; misses are prepared, dispatched
    through the executor, and gathered with per-future exception
    capture so no single failure aborts the batch.  Entries are
    returned in unit order.

    When ``executor`` is None one is chosen from the flight:
    a ``ParslExecutor`` if it carries a ``parsl_config``, else a
    ``LocalExecutor``.  A caller may pass an executor explicitly
    (tests do, to pin the path)."""
    validate_flight(flight)
    os.makedirs(flight.root, exist_ok=True)
    serialize_flight(flight)

    owns_executor = executor is None
    if executor is None:
        if flight.parsl_config is not None:
            executor = ParslExecutor(flight.parsl_config)
        else:
            executor = LocalExecutor()

    results = [None] * len(flight.units)
    try:
        # Pass 1: report hits immediately, dispatch misses.
        pending = []                  # (index, unit, future)
        for index, unit in enumerate(flight.units):
            wingbeat_dir = unit_run_dir(flight, unit)
            if is_cache_hit(unit, wingbeat_dir):
                entry = report_entry_from_status(unit, wingbeat_dir)
                results[index] = entry
                _fire(flight, entry)
            else:
                _prepare_miss(flight, unit, wingbeat_dir)
                future = executor.submit_unit(
                    unit, wingbeat_dir, flight.default_wingbeat
                )
                pending.append((index, unit, future))

        # Pass 2: gather the dispatched units.
        for index, unit, future in pending:
            entry = _collect(flight, unit, future)
            results[index] = entry
            _fire(flight, entry)
    finally:
        if owns_executor:
            executor.close()

    return FlightReport(entries=results)
