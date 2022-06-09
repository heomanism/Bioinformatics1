# 생물정보학 및 실습 1 자유 프로젝트

## 프로젝트 목표
1. **Fig S3C 재현해보기** : transcriptome에서 error가 많이 나오는 부분들 주변 서열을 모아서 문맥을 파악
2. **LIN28A binding motif predictor 만들기** : Fig S3C를 만드는 도중 얻을 수 있는 LIN28A-bound sequence들로 기계학습으로 예측 모델을 만들고 평가합니다. CNN이나 random forest, kNN 같은 걸 쓸 수도 있지만 단순하게 PSSM을 보정해서 쓰는 것도 가능함
3. 위의 내용을 우선은 R을 사용하여 분석해보고, 시간이 허락한다면 Python으로도 해볼 예정

## 프로젝트 관련 논문 정리
- Monoclonal anti-LIN28A antibody인 35L33G를 이용하여 CLIP-seq data(CLIP-35L33G.bam)를 얻음 
- CLIP-seq data는 single-nucleotide resolution에서 protein binding sites를 mapping하는 것을 도움
- most substitutions in CLIP tags는 base 'G'에서 대부분 발견되었고, 나머지 'A','C','U'는 RNA-seq library에서의 error frequency와 비교햇을 때 유사한 것을 확인함
- LIN28A-binding sites를 transcriptome-wide level에서 Shannon's entropy information를 통해서 찾음 (특정 position에서 nucleotide composition의 randomness를 정량화 하는 방법)
- Crosslinking-induced reverese-transcription error score(CRES)가 Shannon's entropy로 결정됨
- mutated nucleotide를 포함하는 sequence를 centered 시켜서 WebLogo로 visualization 시킴 (Fig 2A, Fig S3C)
- hexamer seuqence를 binding site를 centered한 것에서 -2 ~ +3으로 설정함


## 프로젝트 내용 
### ~ 5/27
- 프로젝트 설정
- 프로젝트 관련 논문 내용 정리
- LIN28A의 Shannon entropy 플롯을 그림

### ~ 06/03
- 프로젝트 관련 논문 내용 정리 추가
- Transcriptome-wide level에서 CLIP-35L33G.bam으로부터 CLIP-35L33G.pileup을 만든 후, Shannon's entropy information를 구함
- 기계학습 모델인 KNN, SVM, randomForest 개념 공부를 진행하였고, 모델 학습 방법에 대한 공부를 함
- 현재 계획으로는 10(5) fold Cross-validation 과정으로 모델의 성능을 평가할 예정이지만, 문헌 조사를 통해 마땅한 Replication dataset을 찾는다면 데이터셋을 모델에 적용할 예정

### ~ 06/10
- Chromosome이 정해지지 않은 Scaffold를 filtering하였고, 분석 과정에서는 Mitochondria chromosome을 filering하고 진행할 
- Reverse Strand가 pileup 파일에서 소문자로 인식되는 것을 확인하고, upper로 바꾼 다음, entropy를 다시 구함 (기존 코드에는 대소문자를 나눠서 엔트로피를 구하는 것으로 확인하였기 때문)
-
