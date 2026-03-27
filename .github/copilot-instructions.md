# HPC execution policy

All non-trivial commands must follow this workflow:

1. Never run long or compute-heavy commands on the login node.
2. First ensure execution happens inside a persistent `screen` session on the login node.
3. Before running a task, estimate the minimum required walltime, memory, and CPU from:
   - the type of task groups you have to run (e.g., alignment, variant calling, annotation)
   - input size 
   - previous benchmark files if available
   - try to be conservative in your estimation, and report any assumptions you make.
   - increase the estimate if you are uncertain, but do not overestimate by more than 2x.
   - If it crashes, you can always adjust the estimate and rerun the task.
4. Request a compute node or submit a job using those estimated resources.
5. Run each group of tasks in its own isolated execution context:
   - separate `screen` session, or
   - separate scheduler job

## Hard rules

- Never request a node before creating or attaching a `screen` session.
- Never run heavy commands directly in the login shell.
- If resource estimation is uncertain, choose the smallest safe profile and report the assumption.
- Each task must run in its own isolated session or scheduler allocation.