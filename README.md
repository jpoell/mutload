# mutload
Mutation Load Determination from Duplex Consensus Base Coverage

This repository contains the bioinformatics pipeline to derive a measure of mutation load from sequencing reads that have been compiled from a duplex consensus. This will generally only work if sequencing depth was sufficient to cover both the Watson and the Crick strand of most DNA templates. This is generally achieved either by using a small capture panel (Schmitt et al. PNAS 2012) or by bottlenecking the initial input DNA (Hoang et al. PNAS 2016). The enclosed algorithms rely on coverage from multiple samples to filter out problematic base positions, so reduction of complexity by using a capture panel is the way to go here. To accommodate the user, I have included a workflow that derives duplex consensus reads using the original DCS algorithm as best described in Kennedy et al. Nat Protoc 2014. Although not included here, I would like to mention that alternatives are available as well, such as ConsensusCruncher (Wang et al. NAR 2019). 

Once the consensus reads are derived, all bases are piled up using user-defined regions (most likely the applied capture panel). Mutations are scored as base positions with one or more variant base calls. Variant is defined as any other base than the predominant base in that sample (so not necessarily the genome reference base!). Positions that are frequently altered are filtered out. This most notably includes germline single nucleotide polymorphisms, but also positions that are represented by significantly more samples than expected by chance. The number of mutations are reported, and the cumulative base coverage. Mutation load is reported as the number of mutations per million bases covered. 

The suggested workflow to derive duplex consensus reads requires the following software:
- A command-line interpreter with paste, gzip, zgrep, zcat, and perl
- Python 2.7 and the UnifiedConsensusMaker.py script and its dependencies (https://github.com/Kennedy-Lab-UW/Duplex-Sequencing)
- gatk4 (alternatively Picard) for FastqToSam (4.2.4.0)
- bwa (0.7.17)
- samtools (1.12)

Coverage pileup and rare somatic mutation counting requires:
- bcftools (1.12)
- R (3.6 or 4) and the foreach package

Suggested versions are in parentheses. For most software newer versions are expected to work as well, except for UnifiedConsensusMaker.py (also see Notes below). 

Notes:
The Loeb lab has a newer version of their Duplex Consensus Sequencing Pipeline, but that is not the one used here.
In our experience, bwa aln produced lower background than bwa mem.
