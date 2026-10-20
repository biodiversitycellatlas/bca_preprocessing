## Git Submodules

This repository uses several external pipelines and tools as git submodules. These submodules are not redistributed here; they are referenced directly from their official public repositories.


### Included Submodules

| Path                         | Source Repository                                          | Description |
|-----------------------------|--------------------------------------------------------------|-------------|
| `submodules/10x_saturate`   | https://github.com/bonitavw/10x_saturate.git (branch `CR_UR_support`), fork of [zolotarovgl/10x_saturate](https://github.com/zolotarovgl/10x_saturate) | 10x Genomics saturation and library complexity estimation tool. |
| `submodules/GeneExt`        | https://github.com/zolotarovgl/GeneExt.git                  | Gene model adjustment tool for improved scRNA-seq quantification. |
| `submodules/pavianCore`     | https://github.com/Enthusiasm23/pavianCore.git              | Command-line visualization for kraken reports, extension of [Pavian](https://github.com/fbreitwieser/pavian). |
