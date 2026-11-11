# Environment and compatibility report

Audit performed 2026-11-09 America/New_York (UTC log date may be November 10).
The accessible host is a Linux container, **not the user's Windows/WSL2 machine**.
It exposes approximately 9.7 GiB RAM, 60 GiB free disk and no NVIDIA device or
nvidia-smi. No CUDA model, checkpoint, production MD or cloud GPU was run.
See `reports/session_environment.json` for raw read-only diagnostics.
An isolated Python 3.12 virtual environment was created outside the repository;
the existing system Python and CUDA configuration were not changed.

Run this first inside the user's WSL Ubuntu shell (no dependencies required):

```bash
python3 environment/diagnose.py --out reports/wsl_environment.json
```

The reported Ryzen model, 8 GB GTX 1080, available system RAM, Windows NVIDIA
driver and 150 GB storage remain unverified until that report exists.
The diagnostic output is exclusive-create: use a fresh name to repeat it.

| Workload | Intended location | Current evidence |
|---|---|---|
| RDKit preparation, validation, tabulation, dashboard | Local CPU | CPU code tested on this host; target ligand generation gated |
| Original RFdiffusion small monomer | Test isolated older CUDA stack on GTX 1080, or cloud | Official command attempted here; missing torch; no inference |
| LigandMPNN | Local CPU or compatible GPU environment | CPU reference inference and interface checks passed; GPU untested |
| RFdiffusion3 | Prefer modern cloud GPU after exact environment check | Latest pinned source uses modern PyTorch; Pascal unsupported by common wheels/kernels |
| Independent complex prediction | Prefer modern cloud GPU | Must benchmark chosen model and sequence length before scaling |
| Vina docking | Local CPU | Protocol planned; no docking executed |
| OpenMM explicit-solvent MD | GTX 1080 only after compatible-build test, or cloud | Coordinate analysis tested; MD engine/parameters not validated |

GTX 1080 is Pascal, compute capability **6.1**. It lacks native BF16 and Tensor
Cores. CUDA detection alone says nothing about wheel architectures or attention
kernels. NVIDIA lists Pascal's last toolkit support as CUDA 12.x. PyTorch removed
Pascal from its CUDA 12.8/12.9 binary builds starting with PyTorch 2.8. Thus even a
CUDA 12.x label is insufficient: check the exact wheel's `get_arch_list()` and run
a small allocation/matrix operation plus the actual model kernels. Do not install
the newest PyTorch wheel blindly. Legacy Foundry inference disables some fast
kernels but does not make an incompatible PyTorch build support Pascal.

Do not install a Linux display driver inside WSL or replace an existing CUDA/Conda
installation for this project. Create separate environments for original RFdiffusion,
Foundry, LigandMPNN, and analysis. The original published RFdiffusion YAML is
preserved under `configs/reference_examples/SE3nv.original.yml`; its old Python,
PyTorch and CUDA pins are not interchangeable with current Foundry requirements.
Neither model stack has been installed/validated here. Checkpoint hashes are absent
until downloaded; source commit hashes must not be mistaken for checkpoint hashes.

Sources checked:
- https://developer.nvidia.com/cuda/gpus/legacy
- https://docs.nvidia.com/datacenter/tesla/drivers/latest/cuda-toolkit-driver-and-architecture-matrix.html
- https://dev-discuss.pytorch.org/t/cuda-toolkit-version-and-architecture-support-update-maxwell-and-pascal-architecture-support-removed-in-cuda-12-8-and-12-9-builds/3128
- https://github.com/RosettaCommons/foundry/blob/production/pyproject.toml
- https://github.com/RosettaCommons/foundry/blob/production/models/rfd3/docs/inference_acceleration.md

## Cost and resource policy

No provider credentials, rented instance, paid API call, or automatic cloud
provisioning is included. Default spending cap is zero. Model throughput and cost
are **unmeasured**; quote an actual provider rate only when selecting hardware.
Budget = measured pilot GPU-hours × current quoted hourly rate, plus storage and
egress. A planning example of 2 GPU-hours at an assumed $2/hour is $4, not a quote
or prediction. Run one design first, measure memory/runtime, then estimate the
remaining campaign. Provider-side auto-shutdown and billing limits are required;
a local process timeout cannot enforce a provider billing ceiling.
