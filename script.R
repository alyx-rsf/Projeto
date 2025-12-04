# ===============================================================================
# title: "Arquitetura genética da resistência à murcha de *Fusarium* na população 'IAC-Tybatã x Branquinho'de feijoeiro comum"
# subtitle: "Simulação de análise fenotípica"
# author: "Alyx R. S. de Freitas" 
# ====================================Sumário===================================
# 1. Introdução 
# 2. Simulação 
# 3. Análise descritiva 
# 4. Severidade aos 21 DAI
# 5. ANOVA
# 6. Cálculo da Área Abaixo da Curva (AACPD) 
# 7. Genótipos
# 8. Visualização
# 9.Parentais x RILs
# 1. Introdução=================================
# Problema: Murcha de *Fusarium* é uma das principais doenças do feijoeiro
# Objetivo: Simular dados fenotípicos para análise de resistência
# População: Cruzamento entre IAC-Tybatã (resistente) × Branquinho (suscetível)
# Aplicação: Métodos estatísticos para classificação de genótipos
# 2. Simulação=================================

#Carregar pacotes necessários
library(tidyverse)
library(readxl)
library(lme4)
library(agricolae)

#Carregar id dos genótipos
setwd("C:/Documentos/UNICAMP/R/Projeto")
dados_estrutura <- read_excel("id_sementes.xlsx")
#Configurar dados da simulação
set.seed(123)

dados <- dados_estrutura %>%
  mutate(
    perfil = case_when(
      ID == "Branquinho" ~ "suscetivel",
      ID == "Tybatã" ~ "resistente",
      T ~ sample(c("resistente", "moderado", "suscetivel"), n(), replace = T, prob = c(0.2, 0.6, 0.2))))

# Mostrar estrutura dos dados simulados
dados


# Simular severidade com progressão temporal

dados_com_sev <- dados %>%
  mutate(
    SEV = case_when(
      # DAI 7 - Primeira avaliação (sintomas iniciais)
      DAI == 7 & perfil == "resistente" ~ sample(c(1, 1, 1, 1, 1), n(), replace = TRUE),
      DAI == 7 & perfil == "moderado" ~ sample(c(1, 1, 1, 3, 3), n(), replace = TRUE),
      DAI == 7 & perfil == "suscetivel" ~ sample(c(1, 1, 3, 3, 3, 5), n(), replace = TRUE),
      
      # DAI 14 - Segunda avaliação (progressão moderada)
      DAI == 14 & perfil == "resistente" ~ sample(c(1, 1, 3, 3, 3), n(), replace = TRUE),
      DAI == 14 & perfil == "moderado" ~ sample(c(3, 3, 5, 5, 5), n(), replace = TRUE),
      DAI == 14 & perfil == "suscetivel" ~ sample(c(5, 5, 7, 7, 7), n(), replace = TRUE),
      
      # DAI 21 - Terceira avaliação (máxima expressão)
      DAI == 21 & perfil == "resistente" ~ sample(c(1, 3, 3), n(), replace = TRUE),
      DAI == 21 & perfil == "moderado" ~ sample(c(5, 5, 7), n(), replace = TRUE),
      DAI == 21 & perfil == "suscetivel" ~ sample(c(7, 9, 9), n(), replace = TRUE),
      
      TRUE ~ NA_real_
    )
  )


# 3. Análise descritiva=================================
# filtrar dados 21 DAI
sev_21 <- dados_com_sev %>% 
  filter(DAI == 21) %>% 
  mutate(
    ID = factor(ID),
    REP = factor(REP)
  )

# Média de severidade por genotipo
class_sev <- sev_21 %>% 
  group_by(ID) %>% 
  summarise(
    media_SEV = mean(SEV, na.rm = T),
    dp_SEV = sd(SEV, na.rm = T),
    n_rep = n()
  ) %>% 
  mutate(
    Perfil = case_when(
      media_SEV <= 3  ~ "Resistente",
      media_SEV <= 5 ~ "Moderado",
      T ~ "Suscetível"
    ),
    Perfil = factor(Perfil,
                    levels = c("Resistente", "Moderado","Suscetível"))
  )
# Mostrar tabela de classificação
head(class_sev, 10)

# 4. Visualização da severidade aos 21 DAI=================================
ggplot(data = sev_21, aes(x = "", y = SEV)) +
  geom_boxplot(fill = "lightgreen", color = "darkgreen") +  
  geom_jitter(aes(color = ID), width = 0.1, alpha = 0.3, size = 2) +
  labs(
    title = "Severidade Fenotípica aos 21 DAI",
    subtitle = "Distribuição por Genótipo",
    x = "Todos os Genótipos",
    y = "Severidade (0-9)") +
  theme_minimal() +
  theme(legend.position = "none")
# 5. ANOVA =================================
# ANOVA para severidade


anova_sev <- aov(SEV ~ ID, data = sev_21)
anova_sev_sum <- summary(anova_sev)
summary(anova_sev) #??????????????????

shapiro.test(residuals(anova_sev))

# Teste de Tukey 
tuckey_sev <- HSD.test(anova_sev, trt = "ID")
tuckey_sev$groups

# 6. Cálculo da Área Abaixo da Curva (AACPD)=================================
# AACPD é uma medida que integra a severidade da doença ao longo do tempo.
#Quanto maior a AACPD, mais severa a progressão.

#AACPD não funcionou com o data frame completo:
dados_temp <- dados_com_sev %>%
  select(ID, REP, DAI, SEV)

aacpd <- dados_temp %>%
  group_by(ID, REP) %>%
  summarise(
    AACPD = audpc(evaluation = SEV, dates = DAI),
    .groups = "drop"
  )

# Estatísticas AACPD
estatisticas_aacpd_geral <- aacpd %>%
  summarise(
    n_plantas = n(),
    n_genotipos = n_distinct(ID),
    media_AACPD = mean(AACPD, na.rm = T),
    mediana_AACPD = median(AACPD, na.rm = T),
    dp_AACPD = sd(AACPD, na.rm = T),
    min_AACPD = min(AACPD, na.rm = T),
    max_AACPD = max(AACPD, na.rm = T),
    q1 = quantile(AACPD, 0.25, na.rm = T),
    q3 = quantile(AACPD, 0.75, na.rm = T)
  )
print(estatisticas_aacpd_geral)

ggplot(data = aacpd_clean, aes(x = "", y = AACPD)) +
  geom_boxplot(fill = "darkorchid", color = "darkorchid4") +  
  geom_jitter(aes(color = ID), width = 0.1, alpha = 0.3, size = 2) +
  labs(
    title = "AACPD por genótipo",
    x = "Todos os Genótipos",
    y = "AACPD (0-100)") +
  theme_minimal() +
  theme(legend.position = "none")

#ANOVA?
aacpd_clean <-  aacpd %>%
  filter(!is.na(AACPD))

anova_resultado <- aov(AACPD ~ ID, data = aacpd_clean)
anova_summary <- summary(anova_resultado)
anova_summary

tuckey_sev <- HSD.test(anova_resultado, trt = "ID")
tuckey_sev$groups

# Coeficiente de variação experimental (com arredondamento)
cv_experimental <- round((sqrt(anova_summary[[1]]$`Mean Sq`[2]) / mean(aacpd_clean$AACPD)) * 100, 1)
cv_experimental

teste <- compare_means(AACPD ~ ID,  data = aacpd_clean)
teste

# 7. Genótipos =================================
estatisticas_genotipos <- aacpd %>%
  group_by(ID) %>%
  summarise(
    n_rep = n(),
    media_AACPD = mean(AACPD, na.rm = TRUE),
    dp_AACPD = sd(AACPD, na.rm = TRUE),
    cv_AACPD = (sd(AACPD, na.rm = TRUE) / mean(AACPD, na.rm = TRUE)) * 100,
    .groups = "drop"
  )
head(estatisticas_genotipos, 10)

# Classificação
classificacao <- estatisticas_genotipos %>%
  mutate(
    Classificacao = case_when(
      media_AACPD < 30 ~ "Altamente Resistente",
      media_AACPD < 50 ~ "Resistente", 
      media_AACPD < 70 ~ "Moderado",
      media_AACPD < 90 ~ "Suscetível",
      TRUE ~ "Altamente Suscetível"
    ),
    Classificacao = factor(Classificacao, 
                           levels = c("Altamente Resistente", "Resistente", 
                                      "Moderado", "Suscetível", 
                                      "Altamente Suscetível"))
  )

# 8. Visalização =================================

# Histograma da distribuição da AACPD

aacpd_clean <-  aacpd %>%
  filter(!is.na(AACPD))

ggplot(aacpd_clean, aes(x = AACPD)) +
  geom_histogram(bins = 12, fill = "steelblue", color = "white", alpha = 0.7) +
  geom_vline(xintercept = mean(aacpd_clean$AACPD), 
             color = "red", linetype = "dashed", size = 1) +
  geom_vline(xintercept = median(aacpd_clean$AACPD), 
             color = "blue", linetype = "dotted", size = 1) +
  labs(title = "Distribuição da AACPD",
       subtitle = paste("Média =", round(mean(aacpd_clean$AACPD), 1), 
                        "| Mediana =", round(median(aacpd_clean$AACPD), 1)),
       x = "AACPD",
       y = "Frequência") +
  theme_minimal()

# Boxplot da AACPD por genótipo (ordenado) (ficou horrivel)
ggplot(aacpd_clean, aes(x = reorder(ID, AACPD), y = AACPD)) +
  geom_boxplot(fill = "lightblue", alpha = 0.7) +
  labs(title = "Distribuição da AACPD por Genótipo",
       subtitle = "Ordenado pela mediana da AACPD",
       x = "Genótipo",
       y = "AACPD") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6))

# 9. Parentais x RILs ==========================================

parentais_data <- aacpd_clean %>%
  filter(ID %in% c("Branquinho", "Tybatã"))

ggplot(parentais_data, aes(x = ID, y = AACPD, fill = ID)) +
  geom_boxplot(alpha = 0.7) +
  labs(title = "Comparação da AACPD entre Parentais",
       x = "Parental",
       y = "AACPD") +
  scale_fill_manual(values = c("Branquinho" = "red", "Tybatã" = "green")) +
  theme_minimal()

# RILs

# 10 genótipos mais resistentes (menor AACPD)
top10_resistentes <- classificacao %>%
  arrange(media_AACPD) 

# 10 genótipos mais suscetíveis (maior AACPD)
top10_suscetiveis <- classificacao %>%
  arrange(desc(media_AACPD)) 

top10_resistentes

top10_suscetiveis
