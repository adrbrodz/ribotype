# Ribotype

**Ribotype** was developed as a computational method
for comprehensive annotation of riboswitch regulatory units, including aptamer,
expression platform and prediction of mechanism of riboswitch action. The pipeline implements
a range of bioinformatic tools to:
<ul>
<li>identify hallmarks of expression platform switching, based on analysis of
secondary structure dynamics</li>
<li>identify RBS sequences in bacterial UTR regions</li>
<li>obtain covariance profiles for known riboswitch classes</li>
<li>identify rho-independent transcription termination signals in bacterial genomes</li>
</ul>

## Bioinformatic pipeline
<p align="center">
  <img src="https://github.com/adrbrodz/adrbrodz.github.io/blob/main/src/ribotype-flow2.png?raw=true" width="30%" height="30%">
</p>
<p align = center>
  Flowchart representing the general steps in the identification of </br>
riboswitch sequences and their functional structures.</p>

## Example input

<p align="center">
  <img width="500" height="200" alt="Zrzut ekranu 2025-09-10 102715" src="https://github.com/user-attachments/assets/de8335e0-b31b-4c29-b4cf-db536a264613" />
</p>
<p align="center">
  Complete genomic sequence for <i>Klebsiella pasteurii</i> in the FASTA format,<br/>spanning ~6 million nucleotides.
</p>

## Example output

<p align="center">
  <img width="500" height="375" alt="Zrzut ekranu 2025-09-10 103018" src="https://github.com/user-attachments/assets/2806b117-30a4-45ed-b5b9-1299e5d46e94" />
</p>
<p align="center">
The example output of automatic riboswitch annotation and mechanism<br/>
prediction for <i>Bacillus subtilis</i>.
</p>


## Data visualization

<p align="center">
  <img width="350" height="600" alt="Zrzut ekranu 2025-09-10 101755" src="https://github.com/user-attachments/assets/a2786633-4a35-4e17-a4b3-896a0865f97b" />
</p>
<p align = center>
The visualization of the identified purine riboswitch representing the
mechanism of transcription termination;<br/>(A) — transcription termination,
(B) — transcription proceeded normally.</p>
<br/><br/>
<p align="center">
  <img width="350" height="600" alt="Zrzut ekranu 2025-09-10 102257" src="https://github.com/user-attachments/assets/f9f21ee9-4793-4e83-8ed2-f3fb5d6bedd1" />
</p>
<p align = center>
The visualization of the identified SAM riboswitch representing the
mechanism of translation inhibition;<br/>(A) - translation proceeded normally,
(B) – translation inhibition.
</p>
