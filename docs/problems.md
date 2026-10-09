# Problems log

Record each problem: what happened, why, what you tried, what fixed it.
These become evidence for the Outcome and Reflection chapters.

## 2026-10-09: Synthetic test genome too short for Prodigal training
- **What:** Pyrodigal warned that the 60 kb test sequence was below the recommended 100 kb for training.
- **Fix:** Below 100 kb, CAM switches to metagenomic mode automatically and logs which mode it used.
