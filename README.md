# Preprocess eCLIP data from Munschauer lab

The eCLIP libraries prepared in the Munschauer lab have a 10 base pair UMI in the first bases of R2. This is a simple script that trims the UMI and puts it in the sequence header so that one can align and then perform UMI-aware deduplicating (e.g. with Picard tools). 

Libraries are reported here: https://www.biorxiv.org/content/10.1101/2020.07.15.204404v1


MWE: 
```
python CLIP_trim.py -a data/CLIP_1.fastq.gz -b data/CLIP_2.fastq.gz -o trimmed_UMI
```
<br><br>
