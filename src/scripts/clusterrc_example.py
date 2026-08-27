#!/usr/bin/env python3

"""clusterrc.py -- starter written by cluster_probe.py.

Review every value below.  A FILL IN setting could not be
discovered -- complete it from your site's documentation.
"""


def parameters_and_defaults():
    """Return this site's cluster dispatch settings."""
    return {
        # The scheduler queues you can submit to (first is default).
        'partitions': [
            'rulisp-lab', 'general', 'requeue'
        ],
        # Shell commands that set up Imago before a job runs.  A
        # SLURM worker runs these in a NON-interactive shell, so the
        # full activation chain is spelled out here (the interactive
        # cpg/simago aliases do not exist there).  Order matters:
        # conda.sh enables `conda activate` non-interactively; the cpg
        # env + its LD_LIBRARY_PATH supply the HDF5 runtime the imago
        # Fortran binary links against; the venv supplies parsl / numpy
        # / h5py / ase; imagorc is sourced LAST because the venv's
        # activate overwrites PYTHONPATH, and imagorc must prepend
        # $IMAGO_BIN to it so `import imago` / kaleidoscope resolve.
        #
        # The venv is STAGED onto the node's local disk rather than
        # activated where it lives on the share.  Reading it over the
        # shared file system costs roughly half a second per file when
        # the cluster is busy, and starting Python opens thousands of
        # files, so `import parsl` has been measured at over ten
        # minutes of waiting against one second of CPU.  Parsl allows
        # its interchange helper a hard-coded 120 seconds to start, so
        # a slow share does not merely delay a campaign -- it stops one
        # from starting at all.  stage_venv.sh streams the same
        # environment as a single tar file, which costs seconds, and
        # Parsl finds the staged interchange.py and
        # process_worker_pool.py because it launches both by bare name
        # through PATH.  Should staging fail -- a full scratch disk,
        # say -- the `||` falls back to the shared copy, so a job runs
        # slowly rather than not at all.
        'worker_init': [
            'source /cluster/software/common/mamba/1.4.2/'
            'etc/profile.d/conda.sh',
            'conda activate cpg',
            'export LD_LIBRARY_PATH=/cluster/VAST/rulisp-lab/'
            'cpg/mamba/envs/cpg/lib',
            'CPG_IMAGO_VENV=$(/cluster/VAST/rulisp-lab/cpg/bin/'
            'stage_venv.sh imago) '
            '&& source "$CPG_IMAGO_VENV/bin/activate" '
            '|| source /cluster/VAST/rulisp-lab/cpg/'
            'virtual_envs/imago/bin/activate',
            'source $HOME/imago/.imago/imagorc',
            'export OMP_NUM_THREADS=1',
        ],
        # The account charged for compute time, if one is required.
        # You may use: rulisp-lab, general.
        'account': 'rulisp-lab',
        # How many CPU cores one node has.
        # Nodes vary -- core counts seen: 36, 40, 48, 64, 96, 104, 112, 128,
        # 256. Choose one.
        'cores_per_node': 40,  # FILL IN
        # Calculations to run at once on a node (blank = auto).
        'workers_per_node': 20,
        # How many cores each calculation uses (one for now).
        'cores_per_worker': 1,
        # How many nodes to request per job.
        'nodes': 1,
        # The most time one job may run, as HH:MM:SS.
        'walltime': '02:00:00',
        # Default work spread when no --dispatch flag is given.
        #   The valid dispatch shapes are:
        ##    'slurm-pooled'  -- one allocation whose packed workers
        ##                       stream many units; best for many
        ##                       small, similar units (e.g. seeding).
        ##    'slurm-per-job' -- one scheduler submission per unit;
        ##                       best for large or heterogeneous units.
        ##    'local'         -- run in-process, no scheduler; the
        ##                       deliberate opt-out, usually given
        ##                       explicitly as --dispatch local.
        'default_topology': 'slurm-pooled',
        # Most allocations to grow to at once when work piles up.
        'max_blocks': 16,
        # A node's physical memory, in megabytes.  This is a capacity
        # figure (what the hardware has), NOT a per-job request: it is
        # not spent as `--mem`, and is reserved as a ceiling for future
        # packing and memory-estimation checks.
        # Nodes vary -- memory sizes (MB) seen: 178729, 372254, 372257,
        # 502123, 502184, 502216, 502256, 502258, 502700, 2050457,
        # 2050460, 2050462.
        'memory_per_node': 502123,
        # Memory one calculation needs, in gigabytes.  This is the
        # per-job request: it becomes SLURM's `--mem` (times the workers
        # packed on a node under the pooled shape).  A future
        # memory-estimator will refine it per structure.
        'memory_per_worker': 10,
        # How one calculation starts ('single' = a serial run).
        'launcher': 'single',
        # Parallel pieces per calculation (one for now).
        'ranks_per_worker': 1,
        # Threads per parallel piece (one for now).
        'threads_per_rank': 1,
        # How to pin a calculation to particular cores (off for now).
        'binding': None,
        # Fine control over where threads run (off for now).
        'omp_places': None,
        # Fine control over how threads are pinned (off for now).
        'omp_proc_bind': None,
        # How many GPUs per node to request (zero = none).
        # Nodes vary -- GPU counts seen: 0, 1, 2, 3, 4, 8, 28. Choose one.
        'gpus_per_node': 0,  # FILL IN
        # Settings that should differ for one particular queue.
        'queue_overrides': {},
        # Named setting groups, one per cluster you use.
        'profiles': {},
        # Extra raw scheduler directives to add.
        'extra_scheduler_options': [],
        # Resources for a driver job that prepares work and hands it
        # out, as opposed to one that runs a calculation.  Used only
        # by `--submit`, which wraps the producer in its own batch
        # job.  Any key here is overridable for a single run with
        # --orchestrator-cores / --orchestrator-memory /
        # --orchestrator-walltime, so this is a default shape rather
        # than a fixed one.
        #
        # The right size depends on the dispatch shape.  Under
        # 'slurm-per-job' or 'slurm-pooled' the driver only builds
        # inputs and submits jobs, so a core or two and a little
        # memory are ample -- what it needs is a walltime long enough
        # to outlast the whole flight it supervises.  Under
        # '--dispatch local' the driver runs every SCF itself, in its
        # own process and one after another, so it needs what a
        # single calculation needs: at least `memory_per_worker`
        # above (10 GB), and enough walltime for every calculation in
        # the manifest end to end.  The values below cover both, with
        # the memory sized for the local case and the walltime at the
        # 48-hour ceiling most compute nodes here enforce.
        'orchestrator': {
            'cores': 2,
            'memory': '16G',
            'walltime': '48:00:00',
        },
        # Resources for a molecular-dynamics job -- an outside program
        # run as many parallel ranks filling one node.  Used by
        # condense.py, which writes the LAMMPS submission file from
        # these settings rather than from a template of its own.
        #
        # Blank ranks means "use a whole node's cores", so this job
        # takes its width from `cores_per_node` above (40) and follows
        # it if that is ever corrected.  The init lines bring LAMMPS
        # within reach the way `worker_init` does for Imago: the group
        # build lives outside the module tree the cluster provides, so
        # its directory is added first, and the module is named
        # `cpg_lammps` rather than `lammps` because the site already
        # ships ten `lammps/*` modules of its own.  Loading it also
        # pulls in the matching OpenMPI, which is what supplies
        # `mpirun` -- an MPI program must be launched by the same MPI
        # stack it was linked against, so the version is pinned there
        # rather than chosen here.
        'md': {
            'ranks': None,
            'walltime': '01:00:00',
            'memory': None,
            'init': [
                'module use /cluster/VAST/rulisp-lab/cpg/modulefiles',
                'module load cpg_lammps',
            ],
        },
    }

