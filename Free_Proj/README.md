# 생물정보학 및 실습 1 자유 프로젝트

## 프로젝트 목표
1. **Fig S3C 재현해보기** : transcriptome에서 error가 많이 나오는 부분들 주변 서열을 모아서 문맥을 파악
2. **LIN28A binding motif predictor 만들기** : Fig S3C를 만드는 도중 얻을 수 있는 LIN28A-bound sequence들로 기계학습으로 예측 모델을 만들고 평가합니다. CNN이나 random forest, kNN 같은 걸 쓸 수도 있지만 단순하게 PSSM을 보정해서 쓰는 것도 가능함

## 프로젝트 관련 논문 정리
- Monoclonal anti-LIN28A antibody인 35L33G를 이용하여 CLIP-seq data(CLIP-35L33G.bam)를 얻음 
- CLIP-seq data는 single-nucleotide resolution에서 protein binding sites를 mapping하는 것을 도움
- most substitutions in CLIP tags는 base 'G'에서 대부분 발견되었고, 나머지 'A','C','U'는 RNA-seq library에서의 error frequency와 비교햇을 때 유사한 것을 확인함
- LIN28A-binding sites를 transcriptome-wide level에서 Shannon's entropy information를 통해서 찾음 (특정 position에서 nucleotide composition의 randomness를 정량화 하는 방법)
- Crosslinking-induced reverese-transcription error score(CRES)가 Shannon's entropy로 결정됨
- mutated nucleotide를 포함하는 sequence를 centered 시켜서 WebLogo로 visualization 시킴 (Fig 2A, Fig S3C)
- hexamer seuqence를 binding site를 centered한 것에서 -2 ~ +3으로 설정함


## 프로젝트 내용 
### ~ 5/27 - python
- 프로젝트 설정
- 프로젝트 관련 논문 내용 정리
- LIN28A의 protein이 얼마나 자기 자신의 서열에 붙는지 궁금해서, LIN28A의 Shannon entropy 플롯을 그림

### ~ 06/03 - python
- 프로젝트 관련 논문 내용 정리 추가
- Transcriptome-wide level에서 CLIP-35L33G.bam으로부터 CLIP-35L33G.pileup을 만든 후, Shannon's entropy information를 구함
- 기계학습 모델인 KNN, SVM, randomForest 개념 공부를 진행하였고, 모델 학습 방법에 대한 공부를 함
- 현재 계획으로는 10(5) fold Cross-validation 과정으로 모델의 성능을 평가할 예정이지만, 문헌 조사를 통해 마땅한 Replication dataset을 찾는다면 데이터셋을 모델에 적용할 예정

### ~ 06/10
- Chromosome이 정해지지 않은 Scaffold를 filtering하였고, 분석 과정에서는 Mitochondria chromosome을 filering하고 진행함 
- (-) Strand가 pileup 파일에서 소문자로 인식되는 것을 확인하고, (+) strand와 (-) strand를 각각 나눠서, entropy를 구함
- Motif의 range는 -10 ~ 10로 설정하였고,  Entropy > 0.8 & reads counts > 이상인 것을 Motif, Entropy < 0.8 & reads counts > 50 이상인 것을 Non-Motif라고 함
- Motif 서열을 뽑아내어 Weblogo 웹사이트를 통해 bits plot, probability plot을 그림
- prediction에 사용된 motif의 범위는 -2 ~ 3임 
- Motif의 수보다, Non-motif의 수가 훨씬 많기 때문에, modeling을 진행하기 위해서 Motif 수 : Non-motif 수 = 1:1이 되도록, Non-motif에서 sampling을 진행함
- 여러가지 기계학습 방법들 (KNN, RandomForest, SVM_Linear, SVM_Radial, Logistic Regression, gradient boosting model)로 predictor를 만든 후, 각각에 대한 성능을 내부 데이터 셋에서 10-fold cross validation으로 비교함, 이때, train set과 test set은 독립적임 

### 프로젝트 고찰
- (+) , (-) strand를 처리하는데 시간이 굉장히 많이 소요됨, srandness에 대한 완벽한 이해가 필요함
- 유명한 기계학습 방법들을 가지고, 성능을 테스트하였지만, 정작 내부 알고리즘을 자세히 모르기 때문에 깊은 이해가 필요함
- Non-motif와 Motif 수의 차이가 너무 많이 나기 때문에(Unbalanced Data), Sampling을 하는 방법 말고 다른 방법이 적용되면 더 좋겠다는 생각을 함
- Prediction에 사용된 motif의 범위가 -2 ~ 3 인데, 특정 bp 씩 늘려가면서 predictor를 test하는 것이 필요함, 현재 생각으로는 bp를 늘리면 늘릴수록 성능이 안좋아질 것 같음
