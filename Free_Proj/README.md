# 생물정보학 및 실습 1 자유 프로젝트

## 프로젝트 목표
1. **Fig S3C 재현해보기** : transcriptome에서 error가 많이 나오는 부분들 주변 서열을 모아서 문맥을 파악
2. **LIN28A binding motif predictor 만들기** : Fig S3C를 만드는 도중 얻을 수 있는 LIN28A-bound sequence들로 기계학습으로 예측 모델을 만들고 평가합니다. CNN이나 random forest, kNN 같은 걸 쓸 수도 있지만 단순하게 PSSM을 보정해서 쓰는 것도 가능함
3. **2차구조를 같이 고려한 LIN28A binding motif predictor 만들기** : 윗 주제에 2차구조까지 더 가미해서 정확도를 높임
4. 위의 내용을 우선은 R을 사용하여 분석해보고, 시간이 허락한다면 Python으로도 해볼 예정

## 프로젝트 관련 논문 정리
- CLIP-seq data는 single-nucleotide resolution에서 protein binding sites를 mapping하는 것을 도움
- most substitutions in CLIP tags는 base 'G'에서 대부분 발견되었고, 나머지 'A','C','U'는 RNA-seq library에서의 error frequency와 비교햇을 때 유사한 것을 확인함
- LIN28A-binding sites를 transcriptome-wide level에서 Shannon's information entropy를 통해서 찾음 (특정 position에서 nucleotide composition의 randomness를 정량화 하는 방법)


## 프로젝트 내용 
1번 2번 3번의 목표를 순차적으로 달성할 예정
### ~ 5/27
- 프로젝트 설정
- 프로젝트 관련 논문 내용 정리
- LIN28A의 Shannon entropy 플롯을 그


