#!/bin/bash
if find logfile_non-ref.txt;
then
    cat logfile_non-ref.txt >> logfile_non-ref_prev.txt
fi
${HOME}/github-repos/scripts/pipeline_non-ref_tree_build_1.4-dockerized.sh |& tee logfile_non-ref.txt
