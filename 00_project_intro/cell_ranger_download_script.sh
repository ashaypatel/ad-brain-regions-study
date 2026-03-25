#!/bin/bash


# Make sure you are in your home directory on TSCC (or desired directory)
# For more info: https://www.10xgenomics.com/support/software/cell-ranger/latest
# For more info: https://www.10xgenomics.com/support/software/cell-ranger-atac/downloads

# For single cell/nuclear RNAseq
curl -o cellranger-9.0.1.tar.xz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-9.0.1.tar.xz?Expires=1758109891&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=feJZQ8uJRFlRkQzxDqi6-Xl6tzjfG6OfKDQmVy8Le06VRHjni96T2E9~B9n8HiuRJ7aZ6JTbctMucPB9V-W6hD9hrhkIvQ5fQTrtxSn-oFGJofvPrk7F66lf650nCY266OjCNMc5rqesR2UdmYAtaVuFy07QA3bw0yOCs5MtJqRRXDW~it01WpxaELmTJgYBx2bbX4Ln7DWgXhiOyRxSPRfEvkAGUWvMsEO3aMZ07jkF-XwUF1NYI28RPsPZaIxd3VifeHeWUy9AamMK0ytgjHkcGLCK4oSLalpyJ4OChCMRw5KK-BbGoTqO70xJhcO4cR92KwNpPQmNprFm1vGCcw__"

tar -xf cellranger-9.0.1.tar.xz

curl -O "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz"

tar -xf refdata-gex-GRCh38-2024-A.tar.gz

rm cellranger-9.0.1.tar.xz refdata-gex-GRCh38-2024-A.tar.gz






# For single cell/nuclear ATAC seq
curl -o cellranger-atac-2.2.0.tar.gz "https://cf.10xgenomics.com/releases/cell-atac/cellranger-atac-2.2.0.tar.gz?Expires=1758103056&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=jFOciwmi5A17HRREL3jSBFAtFcFYCiS0DkHOYrFn7T4~k9qq6Y28Hleplbq2p9Xp4iU1wYiS-jOZdjwlgWXogULeuO8hSxPgtHM0Wds-Div6yT8HQ7aOejSy4HcFisGjuq1Yj5GLLmXII~Ugyg4rpdHzpWaMKq6V4sXBP2GLmG5bW2mUnOKJ-Vt-FfuaZl2x9zVrvs87SquAKu8TbYtZGlvt~iQgNe6Hah9pXdjcSnD8LHxyfL~tqcI8023Dtumulovzid3GTWtBKUOvrDeLJVhfqxzkvQ48ud~1ulAr~xRZ0QmAmO6OvaiXsQcWsrv7qkMPgObFXvw2lu88FLYsMQ__"

tar -xf cellranger-atac-2.2.0.tar.gz

curl -O "https://cf.10xgenomics.com/supp/cell-arc/refdata-cellranger-arc-GRCh38-2024-A.tar.gz"

tar -xf refdata-cellranger-arc-GRCh38-2024-A.tar.gz

rm cellranger-atac-2.2.0.tar.gz refdata-cellranger-arc-GRCh38-2024-A.tar.gz

