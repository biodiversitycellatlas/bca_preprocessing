## Git Submodules

This repository uses several external pipelines and tools as git submodules. These submodules are not redistributed here; they are referenced directly from their official public repositories.


### Included Submodules

| Path                         | Source Repository                                          | Description |
|-----------------------------|--------------------------------------------------------------|-------------|
| `submodules/10x_saturate`   | https://github.com/zolotarovgl/10x_saturate.git             | 10x Genomics saturation and library complexity estimation tool. |
| `submodules/GeneExt`        | https://github.com/zolotarovgl/GeneExt.git                  | Gene extension utilities used for reference preprocessing. |
| `submodules/pavianCore`     | https://github.com/Enthusiasm23/pavianCore.git              | Command-line visualization for kraken reports, extension of Pavian. |

---

### Initializing Submodules

After cloning this repository, retrieve all external dependencies with:

```bash
git submodule update --init --recursive
```

To update the submodules to their latest remote commits:
```bash
git submodule update --remote --recursive
```
