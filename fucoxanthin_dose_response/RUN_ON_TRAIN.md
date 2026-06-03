# 기차 안에서 long-run 돌리는 방법

이 v4 버전은 `--save-every` 옵션을 추가했습니다. 계산은 모든 integration step에서 수행하지만, CSV 저장은 일정 간격마다만 합니다.

## 1. 컴파일

```bash
make clean
make
```

## 2. 아주 짧은 long-run 테스트

```bash
make run_long_test
make figures
```

## 3. 기차 안에서 추천하는 moderate run

Mac에서는 노트북이 잠들지 않게 아래처럼 실행하는 것을 추천합니다.

```bash
caffeinate -i bash scripts/run_train_jobs.sh moderate
```

`caffeinate` 없이 실행하려면:

```bash
bash scripts/run_train_jobs.sh moderate
```

이 모드는 다음을 실행합니다.

- `combined_long_mc`
- `fucox0_long`
- `fucox25_long`
- `fucox50_long`
- `fucox100_long`
- figures 자동 생성

기본 설정은 다음입니다.

- `LONG_MC = 200`
- `LONG_STEPS = 100000`
- `LONG_DT = 0.001`
- `LONG_SAVE = 100`

## 4. 더 무거운 full run

```bash
caffeinate -i bash scripts/run_train_jobs.sh full
```

기본 설정은 다음입니다.

- `LONG_MC = 500`
- `LONG_STEPS = 200000`
- `LONG_DT = 0.001`
- `LONG_SAVE = 100`

## 5. 직접 조정해서 실행

```bash
make run_long_mc LONG_MC=500 LONG_STEPS=200000 LONG_DT=0.001 LONG_SAVE=100 LONG_PROGRESS=25
```

```bash
make run_long_fucox LONG_MC=500 LONG_STEPS=200000 LONG_DT=0.001 LONG_SAVE=100 LONG_PROGRESS=25
```

## 6. 진행 상황 확인

로그는 `logs/` 폴더에 저장됩니다.

```bash
tail -f logs/combined_long_mc.log
```

## 7. 결과 확인

주요 결과:

```text
results/combined_long_mc_summary.csv
results/fucox0_long_summary.csv
results/fucox25_long_summary.csv
results/fucox50_long_summary.csv
results/fucox100_long_summary.csv
figure/dose_response/fucoxanthin_dose_response_AUC.png
figure/dose_response/fucoxanthin_dose_response_final_activity.png
```
