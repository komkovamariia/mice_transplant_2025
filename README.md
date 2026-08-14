# TCR-репертуар при аллотрансплантации: воспроизведение и расширение методами вложения последовательностей

Реанализ репертуара α-цепи Т-клеточного рецептора в мышиной модели аллогенной трансплантации костного мозга (BALB/c → C57BL/6). Проект проверяет ранее опубликованные донор-реактивные V-сегменты тремя методологически независимыми линиями анализа и добавляет разрешение на уровне мотивов CDR3 и репертуарной геометрии с помощью библиотеки [mirpy](https://github.com/antigenomics/mirpy).

> *TCR repertoire in an allogeneic bone-marrow transplantation model (BALB/c → C57BL/6): reproduction of the published set-operation / count-based analysis and its extension with sequence-embedding methods (density enrichment, convergent CDR3 motifs, MMD repertoire geometry) from the mirpy toolkit.*

---

## Основной результат

Три независимых подхода — операции над множествами клонотипов, дифференциальный счётный анализ и плотностное вложение последовательностей — **согласованно подтверждают** опубликованные приоритетные V-сегменты:

- **TRAV4-2** — единственный кандидат, занимающий верхние позиции во всех трёх линиях: первое место по плотностному анализу, 4-й ранг по счётному анализу (edgeR) и заметное представительство в witness-анализе;
- **4 из 5** кандидатов статьи (TRAV4-2, TRAV7D-3, TRAV7D-5, TRAV9-4) подтверждены обеими линиями mirpy;
- **TRAV13-1** переопределён как признак аллоответа (обеднён в g1 по обилию) — расхождение содержательно объяснимо и уточняет биологию, а не опровергает её;
- **TRAV6-6**, отнесённый исходной работой к неоднородному фону, на уровне мотивов оказывается наиболее конвергентным (`CALGDMATGGNNKLTF`);
- **устойчивость:** включение J-сегмента не меняет ранжирование V-сегментов (ρ Спирмена = 0,93; лишь 0,11% клонотипов расщепляются по J).

По репертуарной геометрии (MMD до интактного базиса g1): сингенный перенос почти не смещает репертуар (0,069, сопоставимо с внутригрупповой изменчивостью 0,074), тогда как аллогенный перенос смещает его на межлинейное расстояние (g5: 0,134; g6: 0,132).

---

## Структура репозитория

```text
├── venn_original.ipynb        # Исходный анализ: Venn + edgeR/DESeq2/Fisher (TRA/TRB)
├── mirpy_analysis.ipynb       # Вложение, плотность, мотивы, MMD и robustness
├── results_summary.ipynb      # Сводное сравнение подходов
├── study2_alloreactivity/     # Независимый анализ аллореактивных клонотипов
├── методология_mirpy.docx     # Методология и интерпретация
├── results_tables.xlsx        # Агрегированные таблицы
├── figures/                   # Публикационные рисунки
├── requirements.txt           # Python-only зависимости
└── environment.yml            # Каноническое полное окружение статьи
```

### Тетради

| Тетрадь | Содержание |
|---|---|
| **`venn_original.ipynb`** | Первичный анализ, положенный в основу статьи: клонотип = `(CDR3, V)`, диаграммы Венна пересечений групп, локализация V-сегментов в области «только g1», счётный дифференциальный анализ (edgeR, DESeq2, точный тест Фишера), консенсусные клонотипы. |
| **`mirpy_analysis.ipynb`** | Независимая проверка и расширение: единый базис вложения по объединённому пулу, плотностный анализ обогащения (g1 vs g5+g6), конвергентные мотивы CDR3, репертуарные отпечатки Φ(S) и матрица MMD с PERMANOVA, witness-анализ, расширенные биологические контрасты, проверка устойчивости на клонотипах aaVJ, согласование трёх линий. |
| **`results_summary.ipynb`** | Высокоуровневая сводка: читает готовые таблицы результатов и воспроизводит ключевые выводы без повторного тяжёлого вычисления вложения. |

---

## Группы

| Группа | Описание |
|---|---|
| g1 | Интактные реципиенты BALB/c (целевая группа) |
| g2 | Доноры C57BL/6 |
| g3 | Сингенный перенос костного мозга |
| g4 | Контроль кондиционирования |
| g5 | Аллогенный перенос костного мозга |
| g6 | Аллогенный перенос костного мозга + тимус донора |

---

## Воспроизводимое окружение

Для статьи каноническим является **`environment.yml`**. Он устанавливает Python-стек, Jupyter, `ipykernel`, R/rpy2, Bioconductor `edgeR`/`DESeq2`, `mirpy-lib` и `repseq`. Локальный checkout `~/soft/repseq`, ручной `PYTHONPATH`, отдельная установка `dill` и ручная регистрация Jupyter kernel не нужны.

`repseq` устанавливается непосредственно из GitHub с фиксированного коммита:

```text
https://github.com/mmjmike/repseq.git@1c464120ac0675608178608af531f90fcba17deb
```

Фиксация SHA принципиальна: она делает источник `repseq` однозначным и не позволяет будущим изменениям ветки `main` незаметно менять результаты анализа.

### Чистая установка

```bash
git clone https://github.com/komkovamariia/mice_transplant_2025.git
cd mice_transplant_2025

conda env create -f environment.yml
conda activate mice-transplant-2025
```

Никакой `python -m ipykernel install --user` после создания окружения выполнять не требуется. Пакет `ipykernel`, установленный в самом окружении, предоставляет локальный kernelspec `python3` в `${CONDA_PREFIX}/share/jupyter/kernels/python3`. Все воспроизводимые команды ниже запускают именно этот kernel.

Для уже созданного окружения после изменения `environment.yml`:

```bash
conda env update -n mice-transplant-2025 -f environment.yml --prune
conda activate mice-transplant-2025
```

### Проверка установки перед расчётами

```bash
python - <<'PY'
import sys
import repseq
from repseq import clone_filter, clonosets, clustering, diffexp, intersections
from repseq import io, logo, mixcr, slurm, stats, vdjtools
import dill
import matplotlib_venn
from adjustText import adjust_text
import rpy2.robjects as ro
from rpy2.robjects.packages import isinstalled
import mirpy

print("python:", sys.executable)
print("repseq:", repseq.__file__)
print("edgeR:", isinstalled("edgeR"))
print("DESeq2:", isinstalled("DESeq2"))
assert hasattr(repseq, "__path__"), "repseq must be a package, not src/repseq.py"
assert isinstalled("edgeR")
assert isinstalled("DESeq2")
print("Environment OK")
PY
```

Нормальный путь `repseq` должен вести в установленный package внутри окружения, а не в локальный файл `src/repseq.py`.

Проверить встроенный kernel можно без регистрации в пользовательском Jupyter:

```bash
jupyter kernelspec list
```

При активном `mice-transplant-2025` среди доступных kernelspec должен быть `python3`, находящийся внутри текущего conda environment.

### Python-only установка

`requirements.txt` оставлен для задач, которым не нужны R/Bioconductor-блоки. Для полного запуска `venn_original.ipynb` используйте `environment.yml`.

```bash
python -m pip install -r requirements.txt
```

### Данные

Поместите исходные счётные матрицы и подготовленный `clean_clonotypes_aaV.parquet` в `data/` либо задайте путь через переменную окружения:

```bash
export MIRPY_DATA_DIR=/path/to/data
```

### Порядок запуска

1. `venn_original.ipynb`
2. `mirpy_analysis.ipynb`
3. `results_summary.ipynb`

Для серверного воспроизводимого запуска первой тетради после `conda activate mice-transplant-2025`:

```bash
mkdir -p audit_runs logs
jupyter nbconvert \
  --to notebook \
  --execute venn_original.ipynb \
  --ExecutePreprocessor.kernel_name=python3 \
  --ExecutePreprocessor.timeout=-1 \
  --output-dir audit_runs \
  --output 01_venn_original.executed.ipynb \
  2>&1 | tee logs/01_venn_original.log
```

Альтернатива, которая вообще не зависит от того, какое окружение активно в текущем shell:

```bash
conda run -n mice-transplant-2025 jupyter nbconvert \
  --to notebook \
  --execute venn_original.ipynb \
  --ExecutePreprocessor.kernel_name=python3 \
  --ExecutePreprocessor.timeout=-1 \
  --output-dir audit_runs \
  --output 01_venn_original.executed.ipynb
```

Важно: старый `kernelspec` внутри metadata исходного notebook может содержать историческое имя окружения. Для headless-воспроизведения это не используется: параметр `--ExecutePreprocessor.kernel_name=python3` явно выбирает kernel текущего канонического окружения и не требует пользовательской регистрации kernelspec.

`mirpy_analysis.ipynb` содержит наиболее тяжёлый этап; для него требуется существенно больше оперативной памяти, чем для сводной тетради.

---

## Почему Matplotlib закреплён на 3.10.1

Исходный `venn_original.ipynb` использует исторический вызов `plt.boxplot(..., labels=...)`. В Matplotlib 3.11 аргумент `labels` удалён в пользу `tick_labels`, поэтому без фиксации версии старый исходный notebook падает с `TypeError`. Для воспроизведения исходного кода окружение статьи закрепляет Matplotlib 3.10.1.

---

## Доступность данных

Первичные данные секвенирования получены в рамках исходного исследования. Подготовленные клонотипические таблицы и промежуточные артефакты, необходимые для воспроизведения, предоставляются по запросу / размещаются согласно требованиям исходной публикации. Настоящий репозиторий содержит код, агрегированные таблицы результатов и рисунки.

## Ключевые версии

Python 3.12.13, mirpy-lib 3.4.0, numpy 2.5.1, pandas 3.0.3, polars 1.43.0, scipy 1.18.0, scikit-learn 1.9.0, Matplotlib 3.10.1; `repseq` закреплён на Git commit `1c464120ac0675608178608af531f90fcba17deb`. Полный набор зависимостей указан в `environment.yml`.

---

## Исследование 2 — Поиск аллореактивных клонов (mirpy)

Второй, независимый анализ того же набора данных, целиком на инструментарии mirpy / repseq. Многосигнальный конвейер поиска аллореактивных клонотипов α-цепи TCR: профили разнообразия (числа Хилла), биофизическая сигнатура CDR3, сеть сходства и конвергентные кластеры, вероятность генерации (OLGA), публичные клоны и структура компартментов, аннотация VDJdb. Независимо воспроизводит ранжирование V-сегментов Исследования 1 (ρ Спирмена = 1,00) и расширяет его до уровня отдельных клонов.

См. каталог [`study2_alloreactivity/`](study2_alloreactivity/) — ноутбук, рукопись, таблицы результатов и рисунки.
