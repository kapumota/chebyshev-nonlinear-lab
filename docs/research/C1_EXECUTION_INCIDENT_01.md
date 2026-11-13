### C1 execution incident 01

C1 was started from the frozen preregistration state
`c1-preregistration-2026-11-09-r1`.

Execution stopped at command:

`1b5c5db1bd8933d384ec9c0e9f0d6e35a9fe45059f299017b66225a9f5f052f0`

The PETSc process completed normally and produced a durable receipt with
return code zero. The orchestration validator subsequently rejected the
effective PETSc metadata because it compared the complete string returned
by `PetscGetVersion()` against the frozen semantic version `3.26.0`.

Effective PETSc string:

`PETSc Release Version 3.26.0, unknown`

The correction changes only the operational metadata validator. It extracts
the MAJOR.MINOR.SUBMINOR component and compares it against the already frozen
PETSc version.

No scientific command, solver executable, source commit, parameter, start,
seed, repetition, raw CSV, success/failure result, tuning choice, hypothesis,
or statistical protocol is changed.

Previously completed units are preserved. The interrupted unit is recovered
from its existing raw output and successful durable receipt without selective
rerun.
