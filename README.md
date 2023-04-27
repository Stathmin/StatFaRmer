# Wheat_Nitrates_2022
![Picture of wheat with a subscription "I took a picture of wheat, but it came out grainy" ](https://user-images.githubusercontent.com/55657873/229556499-a82053e9-7247-4642-9f5d-9ccd11e9f611.png)


## Задачи:
воспроизвести методы из [этой статьи](https://www.sciencedirect.com/science/article/pii/S0168945223000730) на [данных цифрового фенотипирования Phenospex](data/project_NO3/2022-03-24-Wheat_NO3_#1(b3-6)_20220426_data.zip) с учётом представленных [групп сравнения](data/project_NO3/groups.xlsx).

## Установка:
```{bash}
git clone --branch daniel https://github.com/Stathmin/primate_metabolomics
cd primate_metabolomics
Rscript main.R
```

## Выбор проекта:
В `main.r` замените название эксперимента на название другой папки эксперимента из `data/`:
```{r}
project <- "project_NO3"
```

## Пример:
![красивый график shiny](src/screenshot.png)
