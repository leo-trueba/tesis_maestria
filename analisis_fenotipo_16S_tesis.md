# Análisis 16S raw

## Librerías necesarias

```{r}

library(phyloseq)
library(tidyverse)
library(paletteer)
library(vegan)
library(patchwork)
library(UpSetR)
library(broom)
library(svglite)
library(DESeq2)
library(ape)
library(ggrepel)
library(ggsignif)
library(sf)

```

## Cargar datos

```{r}

otus <- read_table("tabla_otus_tomate_leo.otu.tsv.bz2")
tax <- read_delim("tax_tomate_leo.csv.bz2", delim = ";")

# Para remplazar los nombres de los géneros, con los Incertae Sedis 
# marcados como la última_categoria_disponible__bacterium

tax_sub <- tax %>%
  mutate(tax_completa = paste(Kingdom, Phylum,
                              Class, Order,
                              Family, Genus, sep = ""),
         genero_sub = ifelse(str_detect(tax_completa, "Incertae_Sedis") == TRUE,
                             str_replace(tax_completa, 
                                         "(^.*?)(?=[kpcofg]__Incertae_Sedis).*$", "\\1") %>% 
                             str_replace(".*__(.*$)", "\\1_bacterium"),
                           str_replace(tax_completa, ".*g__", ""))) %>%
  select(OTU, tax_completa, genero_sub)


info <- read_csv("metadatos_completos.csv") %>% 
  mutate(tratamiento = factor(tratamiento, 
                              levels = c("S2", "S3", "S4", "C75"),
                              ordered = TRUE),
         grupo_replica = factor(grupo_replica, ordered  = TRUE), 
         grupo_generacion = factor(grupo_generacion, ordered = TRUE),
         generacion = factor(generacion, ordered = TRUE,
                             levels = c("S", "G1", "G2", "G3")))

orden_grupos <- info %>%select(generacion, tratamiento, replica) %>%
  arrange(generacion, tratamiento) %>% 
  mutate(factor_grupos = paste(generacion, tratamiento, replica, 
                               sep = "_")) %>% pull(factor_grupos)



arbol <- read.tree("tomate_raw.nwk.bz2")

tomate_raw <- phyloseq(column_to_rownames(otus,var = "cluster")  %>% 
                       otu_table(taxa_are_rows = TRUE),
                     column_to_rownames(tax, var = "OTU") %>% as.matrix() %>% tax_table(), 
                     column_to_rownames(info, var = "muestra") %>% sample_data(), 
                     arbol)

tomate_raw

muestras <- colnames(otus)[-1]

fenotipo <- read_csv("fenotipo_plantas.csv") %>%
  mutate(dias = ifelse(generacion %in% c("G1", "G3"), 60, 
                       ifelse(generacion == "G2", 90, NA)), 
         peso_seco_total = peso_seco_aereo + peso_seco_raiz,
         peso_seco_total_rel = peso_seco_total / dias,
         peso_seco_aereo_rel = peso_seco_aereo / dias, 
         peso_seco_raiz_rel = peso_seco_raiz / dias, 
         diametro_tallo_rel = diametro_tallo / dias, 
         longitud_tallo_rel = longitud_tallo / dias) %>%
  filter(tratamiento != "C100")

```


## Requisitos varios


Paletas de colores para los tratamientos y generaciones

```{r}
  
paleta <- c("S2" = "#33D3D3", "S4" = "#CC10CC", "S3" = "#3D55AC",
            "C75" = "#FF8F2E")

paleta_gens <- c("S" = "#2A363BFF", "G1" = "#019875FF", 
                 "G2" = "#EC5CB5FF", "G3" = "#C0392BFF" )


```

Función para combinar tablas de otus, taxonomía y metadatos. Toma el objeto 
de phyloseq que resulta de una operación como `tax_glom()` o 
`transform_sample_counts()` y lo transforma en una tabla para graficar la 
abundancia relativa. Además del objeto de phyloseq, requiere el nombre de la 
variable en dónde colocar las abundancias. Por defecto usa los metadatos 
y nombres de las muestras como se cargaron en  **Cargar datos**.

```{r}

join_phyloseq <- function(obj_phyloseq, var, metadatos=info, seleccion=muestras) {
  otu_temp <- otu_table(obj_phyloseq) %>% 
    as.data.frame() %>%
    rownames_to_column("OTU")

  tax_temp <- tax_table(obj_phyloseq) %>%
    as.data.frame() %>% 
    rownames_to_column("OTU")

  test <- left_join(otu_temp, tax_temp, by = "OTU") %>% 
    pivot_longer(all_of(seleccion), values_to = var,
                 names_to = "muestra") %>% 
  left_join(metadatos, by = "muestra")
}

```

Reemplazo de los nombres de los tratamientos

```{r}


dict_tratamientos <- c("Control", "S2", "S3", "S4")
names(dict_tratamientos) <- c("C75", "S2", "S3", "S4")


```

## Diversidad alfa

Calcular diversidad alfa, añadiendo el índice de Pielou de acuerdo con 
el manual de vegan

```{r}

diversidad <- estimate_richness(tomate_raw, split = TRUE) %>% 
  rownames_to_column(var = "muestra") %>% tibble() %>% 
  left_join(y = info, by = "muestra") %>%
  mutate(pielou = Shannon / log(Observed))

```

### Pruebas estadísticas

```{r}

medidas_diversidad <- colnames(diversidad)[c(2, 3, 5, 7:10,23)]
sub_pruebas <- diversidad %>%
  filter(generacion != "S") 

anova_diversidad <- list()
for (medida in medidas_diversidad){
  lista_medida <- list()
  anova_medida <- aov(formula = get(medida) ~ generacion * tratamiento, 
                      data = sub_pruebas)

  normalidad <- shapiro.test(anova_medida$residuals)
  varianza <- bartlett.test(formula = get(medida) ~ grupo_generacion, 
                            data = sub_pruebas)

  lista_medida[[1]] <- normalidad
  lista_medida[[2]] <- varianza
  lista_medida[[3]] <- summary(anova_medida)

  valor_p <- (unlist(summary(anova_medida)))[c("Pr(>F)1","Pr(>F)2","Pr(>F)3")]

  if(normalidad$p.value >= 0.05 & varianza$p.value >= 0.05 & sum(valor_p < 0.05) > 0) {
    post_hoc <- TukeyHSD(anova_medida) %>% tidy() %>% filter(adj.p.value < 0.05)
    lista_medida[[4]] <- post_hoc
  }
    anova_diversidad[[medida]] <- lista_medida
}


coords_categorias <- read_csv("coords_categorias.csv") %>% 
	select(coords, grupo)

jerarquias <- data.frame(gen_min = c("S", "S", "S", "G1", "G1", "G2"), 
                         gen_max = c("G1", "G2", "G3", "G2", "G3", "G3"), 
                         jerarquia = seq(1.08, 1.48, .08))


signif_observed <- anova_diversidad[["Observed"]][[4]] %>%
  filter(term == "generacion:tratamiento") %>% 
      separate_wider_delim(cols = contrast, delim = "-", 
                           names = c("comp1", "comp2")) %>% 
      separate_wider_delim(cols = comp1, delim = ":", 
                           names = c("generacion1", "tratamiento1"), cols_remove = FALSE) %>%
      separate_wider_delim(cols = comp2, delim = ":", 
                           names = c("generacion2", "tratamiento2"), cols_remove = FALSE) %>% 
      left_join(coords_categorias, by = join_by("comp1" == "grupo")) %>%
      left_join(coords_categorias, by = join_by("comp2" == "grupo"), suffix = c("1", "2")) %>% 
      mutate(xmin = ifelse(coords1 < coords2, coords1, coords2), 
             xmax = ifelse(coords1 < coords2, coords2, coords1), 
             max = 6000, 
             orden = ifelse(coords1 < coords2, comp1, comp2),
             gen_min = ifelse(coords1 < coords2, generacion1, generacion2), 
             gen_max = ifelse(coords1 < coords2, generacion2, generacion1), 
             tipo = ifelse(tratamiento1 == tratamiento2, "dentro", 
                           ifelse(generacion1 == generacion2, "entre_gens", "invalido")),
             anotacion = ifelse(adj.p.value < 0.05  & adj.p.value >= 0.01, "*", 
                                ifelse(adj.p.value < 0.01 & adj.p.value >= 0.001, "**", "***"))) %>% 
      separate_wider_delim(cols = orden, delim = ":", names = c("orden_gen", "orden_trat"))  %>% 
      mutate(orden_gen = factor(orden_gen, levels = c("S", "G1", "G2", "G3"), ordered = TRUE), 
             orden_trat = factor(orden_trat, levels = c("S2", "S3", "S4"), ordered = TRUE)) %>%
      arrange(orden_gen, orden_trat) %>%
      left_join(jerarquias, by = join_by("gen_min", "gen_max")) %>%
      filter(tipo != "invalido") %>% 
      mutate(y_dentro = max * jerarquia, 
             y_entre = ifelse( tipo == "entre_gens", 
                              seq(from = (unique(max) * 1.5), 
                                  to = unique(max) * (1.5 + (0.15 * sum(tipo == "entre_gens"))),
                                        by = unique(max) * 0.15), NA))


for (medida in medidas_diversidad){
  print(medida)
  lista_medida <- list()
  kruskal_medida <- kruskal.test(formula = get(medida) ~ grupo_generacion, 
                                 data = sub_pruebas)

  lista_medida[[1]] <- kruskal_medida

  if (kruskal_medida$p.value < 0.05) {
    wilcox_par <- pairwise.wilcox.test(x = unlist(sub_pruebas[,medida]), 
                         g = sub_pruebas$grupo_generacion)
    lista_medida[[2]] <- wilcox_par
  }

    kruskal_diversidad[[medida]] <- lista_medida
}


```

### Gráficas diversidad alfa

```{r}


observados <- diversidad  %>% 
  ggplot(aes(x = tratamiento, y = Observed, color = generacion)) +
  geom_boxplot() +
  geom_point(position = position_dodge(width = 0.75)) + 
  scale_color_manual(values = paleta_gens) +
  labs(x = NULL, color = "Generación", title = "A") +
  theme_linedraw() + 
  geom_signif(xmin = signif_observed$xmin, 
              xmax = signif_observed$xmax, 
              y = signif_observed$y_dentro,
              annotation = signif_observed$anotacion, 
              color = "black", textsize = 3) +
  geom_signif(xmin = signif_observed$xmin, 
              xmax = signif_observed$xmax, 
              y = signif_observed$y_entre,
              annotation = signif_observed$anotacion, 
              color = "black", linetype = 2, textsize = 3)

chao <- diversidad  %>% 
  ggplot(aes(x = tratamiento, y = Chao1, color = generacion)) +
  geom_boxplot() +
  geom_point(position = position_dodge(width = 0.75)) + 
  scale_color_manual(values = paleta_gens) +
  labs(x = NULL, color = "Generación", title = "B") +
  theme_linedraw()

shannon <- diversidad  %>% 
  ggplot(aes(x = tratamiento, y = Shannon, color = generacion)) +
  geom_boxplot() +
  geom_point(position = position_dodge(width = 0.75)) + 
  scale_color_manual(values = paleta_gens) +
  labs(x = NULL, color = "Generación", title = "C") +
  theme_linedraw()

simpson <- diversidad  %>% 
  ggplot(aes(x = tratamiento, y = Simpson, color = generacion)) +
  geom_boxplot() +
  geom_point(position = position_dodge(width = 0.75)) + 
  scale_color_manual(values = paleta_gens) +
  labs(x = NULL, color = "Generación", title = "D") +
  theme_linedraw()

observados + chao + shannon + simpson + plot_layout(ncol = 2, guides = "collect")

ggsave("figuras/div_alfa.png", units = "cm", height = 18,
       width = 24, dpi = 300)

```



## Abundancia relativa

### Phylum y Clase

```{r}

abundancia_phylum <- tax_glom(tomate_raw, "Phylum") %>% 
  transform_sample_counts(function(x) x/sum(x))

tabla_phylum <- join_phyloseq(abundancia_phylum, "ab_rel")

abundancia_clase <- tax_glom(tomate_raw, "Class") %>% 
  transform_sample_counts(function(x) x/sum(x))

tabla_clase <- join_phyloseq(abundancia_clase, "ab_rel")



tabla_phylum_clase <- rbind(filter(tabla_clase, Phylum == "p__Pseudomonadota"),
                            filter(tabla_phylum,Phylum != "p__Pseudomonadota")) %>% 
  mutate(categoria = ifelse(Phylum == "p__Pseudomonadota", Class, Phylum)) 

orden_categoria <- tabla_phylum_clase %>%
  group_by(categoria) %>%
  summarise(prom_ab = mean(ab_rel)) %>%
  arrange(desc(prom_ab))

dict_categoria <-str_replace(unique(tabla_phylum_clase$categoria), "[pc]__", "") 
names(dict_categoria) <- unique(tabla_phylum_clase$categoria)

tabla_phylum_clase %>%
  left_join(orden_categoria, by = "categoria") %>%
  ggplot(aes(x = grupo_generacion_muestra,
             y = reorder(categoria, prom_ab),
             fill = ab_rel)) + 
  geom_raster() + 
  scale_fill_gradient(low = "#e4e4e4ff", high = "#052629ff", 
                      na.value = "steelblue4", trans = "log10") +
  facet_grid(~tratamiento, scales = "free_x", space = "free_x") +
  scale_y_discrete(labels = dict_categoria)+
  scale_x_discrete(labels = NULL)+
  labs(y = NULL, x = NULL, fill = "Abundancia\nrelativa") +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = -270, size = 6), 
        axis.text.y = element_text(size = 8), 
        legend.position = "none")

  ggsave("figuras/abundancia_phylum_clase.png", width = 15, 
         height = 20, units = "cm", dpi = 300)

  ggsave("figuras/abundancia_phylum_clase.svg", width = 15, 
         height = 20, units = "cm", dpi = 300)

```


### Abundancias relativas por phylum

Pruebas estadísticas de abundancia por phylum, eligiendo las 10 categorías
más abundantes

```{r}


todos_phyla <- orden_categoria$categoria[1:10]
pruebas_phylum <- list()
for (nivel in todos_phyla) {

  sub_categoria <-  tabla_phylum_clase %>%
    filter(categoria == nivel)

  mod_categoria <- aov(formula = ab_rel ~ generacion * tratamiento, data = sub_categoria)

  signif_categoria <- TukeyHSD(mod_categoria) %>% tidy() %>% 
    filter(adj.p.value < 0.05) %>%
    mutate(categoria = nivel, 
           max = max(sub_categoria$ab_rel))

  pruebas_phylum[[nivel]] <- signif_categoria

}

```

Gráficas de abundancia de  Alfaproteobacteria, Gammaproteobacteria, 
Cyanobaceriota y Planctomycetota

```{r}

phyla_seleccionados <- c("c__Alphaproteobacteria", "c__Gammaproteobacteria", 
                         "p__Cyanobacteriota", "p__Planctomycetota")


plots_phylum_signif <- list()
for (nivel in phyla_seleccionados) {

  dim_sub <- pruebas_phylum[[nivel]] %>%
    filter(term == "generacion:tratamiento") %>% dim()

  if(dim_sub[1] != 0) {

    df_coords_signif <- pruebas_phylum[[nivel]] %>%
      filter(term == "generacion:tratamiento") %>%
      separate_wider_delim(cols = contrast, delim = "-", 
                           names = c("comp1", "comp2")) %>%
      separate_wider_delim(cols = comp1, delim = ":", 
                           names = c("generacion1", "tratamiento1"), cols_remove = FALSE) %>%
      separate_wider_delim(cols = comp2, delim = ":", 
                           names = c("generacion2", "tratamiento2"), cols_remove = FALSE) %>%
      left_join(coords_categorias, by = join_by("comp1" == "grupo")) %>%
      left_join(coords_categorias, by = join_by("comp2" == "grupo"), suffix = c("1", "2")) %>%
      mutate(xmin = ifelse(coords1 < coords2, coords1, coords2), 
             xmax = ifelse(coords1 < coords2, coords2, coords1), 
             orden = ifelse(coords1 < coords2, comp1, comp2),
             gen_min = ifelse(coords1 < coords2, generacion1, generacion2), 
             gen_max = ifelse(coords1 < coords2, generacion2, generacion1), 
             tipo = ifelse(tratamiento1 == tratamiento2, "dentro", 
                           ifelse(generacion1 == generacion2, "entre_gens", "invalido")),
             anotacion = ifelse(adj.p.value < 0.05  & adj.p.value >= 0.01, "*", 
                                ifelse(adj.p.value < 0.01 & adj.p.value >= 0.001, "**", "***"))) %>%
      separate_wider_delim(cols = orden, delim = ":", names = c("orden_gen", "orden_trat")) %>%
      mutate(orden_gen = factor(orden_gen, levels = c("S", "G1", "G2", "G3"), ordered = TRUE), 
             orden_trat = factor(orden_trat, levels = c("S2", "S3", "S4"), ordered = TRUE)) %>%
      arrange(orden_gen, orden_trat) %>%
      left_join(jerarquias, by = join_by("gen_min", "gen_max")) %>%
      filter(tipo != "invalido", gen_min != "S") %>%
      mutate(y_dentro = max * jerarquia, 
             y_entre = ifelse( tipo == "entre_gens", 
                              seq(from = (unique(max) * 1.5), 
                                  to = unique(max) * (1.5 + (0.08 * sum(tipo == "entre_gens"))),
                                        by = unique(max) * 0.08), NA))

      dim_signif <- df_coords_signif %>% dim()

      if (dim_signif[1] != 0) {

        plot_nivel <- tabla_phylum_clase %>%
          filter(categoria == nivel) %>%
          ggplot(aes(x = tratamiento, y = ab_rel, color = generacion)) +
          geom_boxplot() + 
          geom_point(position = position_dodge(width = 0.75)) + 
          scale_color_manual(values = paleta_gens, drop = FALSE) + 
          scale_x_discrete(labels = dict_tratamientos) +
          labs(x = NULL, y = NULL, color = "Generación") + 
          facet_wrap(~ categoria, scales = "free") + 
          theme_linedraw() + 
          geom_signif(xmin = df_coords_signif$xmin, 
                      xmax = df_coords_signif$xmax, 
                      y = df_coords_signif$y_dentro,
                      annotation = df_coords_signif$anotacion, 
                      color = "black", textsize = 3) +
          geom_signif(xmin = df_coords_signif$xmin, 
                      xmax = df_coords_signif$xmax, 
                      y = df_coords_signif$y_entre,
                      annotation = df_coords_signif$anotacion, 
                      color = "black", linetype = 2, textsize = 3)
      } else {

        plot_nivel <- tabla_phylum_clase %>%
          filter(categoria == nivel) %>%
          ggplot(aes(x = tratamiento, y = ab_rel, color = generacion)) +
          geom_boxplot() + 
          geom_point(position = position_dodge(width = 0.75)) + 
          scale_color_manual(values = paleta_gens, drop = FALSE) + 
          scale_x_discrete(labels = dict_tratamientos) +
          labs(x = NULL, y = NULL, color = "Generación") + 
          facet_wrap(~ categoria, scales = "free", ncol = 2) + 
          theme_linedraw() 

      }

  } else {

      plot_nivel <- tabla_phylum_clase %>%
        filter(categoria == nivel) %>%
        ggplot(aes(x = tratamiento, y = ab_rel, color = generacion)) +
        geom_boxplot() + 
        geom_point(position = position_dodge(width = 0.75)) + 
        scale_color_manual(values = paleta_gens, drop = FALSE) + 
        scale_x_discrete(labels = dict_tratamientos) +
        labs(x = NULL, y = NULL, color = "Generación") + 
        facet_wrap(~ categoria, scales = "free", ncol = 2) + 
        theme_linedraw() 
  }

      plots_phylum_signif[[nivel]] <- plot_nivel

}

plot_phyla_seleccionados <- wrap_plots(plots_phylum_signif) + plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "A")

ggsave(plot = plot_phyla_seleccionados, filename = "figuras/abundancia_phyla_seleccionados.png", 
       dpi = 300, units = "cm", height = 20, width = 20)

```


Gráficas para los 10  grupos más abundantes (Suplementaria)

```{r}

plots_phylum_todos <- list()
for (nivel in todos_phyla) {

  dim_sub <- pruebas_phylum[[nivel]] %>%
    filter(term == "generacion:tratamiento") %>% dim()

  if(dim_sub[1] != 0) {

    df_coords_signif <- pruebas_phylum[[nivel]] %>%
      filter(term == "generacion:tratamiento") %>%
      separate_wider_delim(cols = contrast, delim = "-", 
                           names = c("comp1", "comp2")) %>%
      separate_wider_delim(cols = comp1, delim = ":", 
                           names = c("generacion1", "tratamiento1"), cols_remove = FALSE) %>%
      separate_wider_delim(cols = comp2, delim = ":", 
                           names = c("generacion2", "tratamiento2"), cols_remove = FALSE) %>%
      left_join(coords_categorias, by = join_by("comp1" == "grupo")) %>%
      left_join(coords_categorias, by = join_by("comp2" == "grupo"), suffix = c("1", "2")) %>%
      mutate(xmin = ifelse(coords1 < coords2, coords1, coords2), 
             xmax = ifelse(coords1 < coords2, coords2, coords1), 
             orden = ifelse(coords1 < coords2, comp1, comp2),
             gen_min = ifelse(coords1 < coords2, generacion1, generacion2), 
             gen_max = ifelse(coords1 < coords2, generacion2, generacion1), 
             tipo = ifelse(tratamiento1 == tratamiento2, "dentro", 
                           ifelse(generacion1 == generacion2, "entre_gens", "invalido")),
             anotacion = ifelse(adj.p.value < 0.05  & adj.p.value >= 0.01, "*", 
                                ifelse(adj.p.value < 0.01 & adj.p.value >= 0.001, "**", "***"))) %>%
      separate_wider_delim(cols = orden, delim = ":", names = c("orden_gen", "orden_trat")) %>%
      mutate(orden_gen = factor(orden_gen, levels = c("S", "G1", "G2", "G3"), ordered = TRUE), 
             orden_trat = factor(orden_trat, levels = c("S2", "S3", "S4"), ordered = TRUE)) %>%
      arrange(orden_gen, orden_trat) %>%
      left_join(jerarquias, by = join_by("gen_min", "gen_max")) %>%
      filter(tipo != "invalido", gen_min != "S") %>%
      mutate(y_dentro = max * jerarquia, 
             y_entre = ifelse( tipo == "entre_gens", 
                              seq(from = (unique(max) * 1.5), 
                                  to = unique(max) * (1.5 + (0.08 * sum(tipo == "entre_gens"))),
                                        by = unique(max) * 0.08), NA))

      dim_signif <- df_coords_signif %>% dim()

      if (dim_signif[1] != 0) {

        plot_nivel <- tabla_phylum_clase %>%
          filter(categoria == nivel) %>%
          ggplot(aes(x = tratamiento, y = ab_rel, color = generacion)) +
          geom_boxplot() + 
          geom_point(position = position_dodge(width = 0.75)) + 
          scale_color_manual(values = paleta_gens, drop = FALSE) + 
          scale_x_discrete(labels = dict_tratamientos) +
          labs(x = NULL, y = NULL, color = "Generación") + 
          facet_wrap(~ categoria, scales = "free") + 
          theme_linedraw() + 
          geom_signif(xmin = df_coords_signif$xmin, 
                      xmax = df_coords_signif$xmax, 
                      y = df_coords_signif$y_dentro,
                      annotation = df_coords_signif$anotacion, 
                      color = "black", textsize = 3) +
          geom_signif(xmin = df_coords_signif$xmin, 
                      xmax = df_coords_signif$xmax, 
                      y = df_coords_signif$y_entre,
                      annotation = df_coords_signif$anotacion, 
                      color = "black", linetype = 2, textsize = 3)
      } else {

        plot_nivel <- tabla_phylum_clase %>%
          filter(categoria == nivel) %>%
          ggplot(aes(x = tratamiento, y = ab_rel, color = generacion)) +
          geom_boxplot() + 
          geom_point(position = position_dodge(width = 0.75)) + 
          scale_color_manual(values = paleta_gens, drop = FALSE) + 
          scale_x_discrete(labels = dict_tratamientos) +
          labs(x = NULL, y = NULL, color = "Generación") + 
          facet_wrap(~ categoria, scales = "free", ncol = 2) + 
          theme_linedraw() 

      }

  } else {

      plot_nivel <- tabla_phylum_clase %>%
        filter(categoria == nivel) %>%
        ggplot(aes(x = tratamiento, y = ab_rel, color = generacion)) +
        geom_boxplot() + 
        geom_point(position = position_dodge(width = 0.75)) + 
        scale_color_manual(values = paleta_gens, drop = FALSE) + 
        scale_x_discrete(labels = dict_tratamientos) +
        labs(x = NULL, y = NULL, color = "Generación") + 
        facet_wrap(~ categoria, scales = "free", ncol = 2) + 
        theme_linedraw() 
  }

      plots_phylum_todos[[nivel]] <- plot_nivel

}

plot_phyla_todos <- wrap_plots(plots_phylum_todos) + plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "A")

ggsave(plot = plot_phyla_todos, filename = "figuras/abundancia_phyla_todos.svg", 
       dpi = 300, units = "cm", height = 50, width = 20, 
       limitsize =FALSE)

```

### Abundancia relativa a nivel de OTU

```{r}

ab_rel_otus <- transform_sample_counts(tomate_raw, function(x) x/sum(x))

tabla_ab_rel_otus <- join_phyloseq(ab_rel_otus, "ab_rel") %>%
  left_join(tax_sub, by = "OTU")

orden_ab_otus <- tabla_ab_rel_otus %>% 
group_by(OTU) %>% 
summarise(prom_ab = mean(ab_rel)) %>% 
arrange(desc(prom_ab))

```

## Diversidad beta

### PCoA con UniFrac

Calclar distancias UniFrac ponderadas y no ponderadas 


```{r}

dist_UniFrac_unw <- UniFrac(tomate_raw, weighted = FALSE)

dist_UniFrac_w <- UniFrac(tomate_raw, weighted = TRUE, normalized = TRUE)
```

```{r}

ord_UniFrac_unw <- ordinate(tomate_raw, method = "PCoA", 
                            distance = dist_UniFrac_unw)

plot_ordination(tomate_raw, ord_UniFrac_unw, color = "generacion", 
                shape = "grupo_tratamiento") +
  geom_point(size = 3) +
  geom_text_repel(aes(label = tratamiento), direction = "both",
            color = "black", size = 2) +
  labs(color = "Generación", shape = NULL) +
  scale_color_manual(values = paleta_gens) +
  scale_shape_manual(values = c(17,16)) +
  theme_linedraw() 

ggsave("figuras/pcoa_UniFrac_unw.png", height = 12, width = 16, units = "cm", 
       dpi = 300)

ggsave("figuras/pcoa_UniFrac_unw.svg", height = 12, width = 16, units = "cm", 
       dpi = 300)

## Prueba estadística

### Efectos separados

adonis2(formula = dist_UniFrac_unw ~ generacion + tratamiento,
	data = info,
	permutations = 9999, by = "margin") %>% tidy()

### Interacción
adonis2(formula = dist_UniFrac_unw ~ generacion:tratamiento,
	data = info,
	permutations = 9999, by = "margin") %>% tidy()


#grupos <- info %>% select(muestra, generacion) %>% pull(generacion)
#betadisper(dist_UniFrac_unw, group = grupos) %>% anova()

ord_UniFrac_w <- ordinate(tomate_raw, method = "PCoA", 
			    distance = dist_UniFrac_w)

plot_ordination(tomate_raw, ord_UniFrac_w, color = "generacion", 
                shape = "grupo_tratamiento") +
  geom_point(size = 4) +
  geom_text_repel(aes(label = tratamiento), direction = "both", force = 0.5,
            color = "black", size = 2) +
  labs(color = "Generacion", shape = NULL) +
  scale_color_manual(values = paleta_gens) +
  scale_shape_manual(values = c(17,16)) +
  theme_linedraw() 

ggsave("figuras/pcoa_UniFrac_w.png", height = 20, width = 22, units = "cm", 
       dpi = 300)

```




### CAP

#### Calcular constriñendo por fenotipo y pH

```{r}

# Quitar muestras de suelo
sub_cap <- subset_samples(tomate_raw, generacion != "S")
#sample_data(sub_cap) <- column_to_rownames(sub_cap, var = "muestra")

# UniFrac Ponderado
UniFrac_cap_w <- UniFrac(sub_cap, weighted = TRUE, normalized = TRUE)

ord_cap_w <- ordinate(sub_cap, method = "CAP", distance = UniFrac_cap_w,
                      formula = ~prom_diametro_tallo + prom_longitud_tallo + 
                        prom_peso_seco_aereo + prom_peso_seco_raiz + prom_pH)

# UniFrac no ponderado

UniFrac_cap_unw <- UniFrac(sub_cap, weighted = FALSE, normalized = TRUE)

ord_cap_unw <- ordinate(sub_cap, method = "CAP", distance = UniFrac_cap_unw,
                        formula = ~prom_diametro_tallo + prom_longitud_tallo + 
                          prom_peso_seco_aereo + prom_peso_seco_raiz + prom_pH)

```

### Pruebas estadísticas

```{r}

# UniFrac Ponderado

muestras_cap <- data.frame(sample_data(sub_cap))

### Efectos separados
adonis2(formula = UniFrac_cap_w ~ generacion + tratamiento,
	data = muestras_cap,
	permutations = 9999, by = "margin") %>% tidy()

### Interacción
adonis2(formula = UniFrac_cap_w ~ generacion:tratamiento,
	data = muestras_cap,
	permutations = 9999, by = "margin") %>% tidy()

# UniFrac No Ponderado

muestras_cap <- data.frame(sample_data(sub_cap))

### Efectos separados
adonis2(formula = UniFrac_cap_unw ~ generacion + tratamiento,
	data = muestras_cap,
	permutations = 9999, by = "margin") %>% tidy()

### Interacción
adonis2(formula = UniFrac_cap_unw ~ generacion:tratamiento,
	data = muestras_cap,
	permutations = 9999, by = "margin") %>% tidy()

```

####  Graficas CAP

```{r}

# Ponderado

loadings_w <- scores(ord_cap_w, display ="bp") %>%
  as.data.frame() %>%
  rownames_to_column("var_fenotipo") %>% tibble() %>% 
  mutate(nombres = c("Diametro tallo", "Longitud tallo", 
                     "Peso seco aereo", "Peso seco raiz", "pH"))

plot_ordination(tomate_raw, ord_cap_w, color = "generacion", 
                shape = "grupo_tratamiento") +
  geom_point(size = 3)+
  geom_text_repel(aes(label = tratamiento), color = "grey60", size = 2) +
  geom_segment(data = loadings_w, 
               aes(x = 0, y = 0, xend = CAP1, yend = CAP2),
               arrow = arrow(length = unit(0.2, "cm")), 
               color = "grey60", inherit.aes = FALSE) +
  geom_text_repel(data = loadings_w, aes(label = nombres, x = CAP1, y = CAP2),
                  inherit.aes = FALSE, min.segment.length = 10000, size = 3) +
  scale_color_manual(values = paleta_gens) +
  scale_shape_manual(values = c(17, 19)) + 
  labs(shape = NULL, color = "Generación") +
  theme_linedraw()

ggsave("figuras/cap_uniFrac_w.png", units = "cm", 
       height = 12, width = 16, dpi = 300)

# No ponderado

loadings_unw <- scores(ord_cap_unw, display ="bp") %>%
  as.data.frame() %>%
  rownames_to_column("var_fenotipo") %>% tibble() %>% 
  mutate(nombres = c("Diametro tallo", "Longitud tallo", 
                     "Peso seco aereo", "Peso seco raiz", "pH"))

plot_ordination(tomate_raw, ord_cap_unw, color = "generacion", 
                shape = "grupo_tratamiento") +
  geom_point(size = 3)+
  geom_text_repel(aes(label = tratamiento), color = "grey60", size = 2) +
  geom_segment(data = loadings_unw, 
               aes(x = 0, y = 0, xend = CAP1, yend = CAP2),
               arrow = arrow(length = unit(0.2, "cm")), 
               color = "black", inherit.aes = FALSE) +
  geom_text_repel(data = loadings_unw, aes(label = nombres, x = CAP1, y = CAP2),
                  inherit.aes = FALSE, min.segment.length = 10000, size = 3) +
  scale_color_manual(values = paleta_gens) +
  scale_shape_manual(values = c(17, 19)) + 
  labs(shape = NULL, color = "Generación") +
  theme_linedraw()


ggsave("figuras/cap_uniFrac_unw.png", units = "cm", 
       height = 12, width = 16, dpi = 300)

ggsave("figuras/cap_uniFrac_unw.svg", units = "cm", 
       height = 12, width = 16, dpi = 300)

```


### PCA fenotipo plantas

```{r}

gen_glom_obj <- tax_glom(tomate_raw, "Genus")

phylum_glom_obj <- tax_glom(tomate_raw, "Phylum")

nmds_genero <- ordinate(gen_glom_obj, method = "NMDS",distance = "bray")

coord_genero <- scores(nmds_genero, display = "species") %>% data.frame()  %>%
	rownames_to_column(var = "OTU") %>%
	left_join(tax_sub, by = "OTU")


coord_muestras_g <- scores(nmds_genero, display = "sites") %>% data.frame() %>% 
	rownames_to_column(var = "muestra") %>% 
	left_join(info, by = "muestra")


abundancia_sub <- prop_total_generos %>%
  arrange(desc(ab_relativa)) %>%
  filter(ab_relativa >= 0.01) %>%
  select(genero_sub, muestra, ab_relativa) %>% pull(genero_sub) %>% unique()

coord_genero %>%
  filter(genero_sub %in% abundancia_sub)


ggplot() + 
	# geom_point(data = coord_genero, aes(x = NMDS2, NMDS2), 
	#            color = "darkorange", alpha = 1.2, pch = 15) +
	# geom_point(data = filter(coord_genero, genero_sub %in% abundancia_sub),
	#             aes(x = NMDS1, NMDS2), 
	#             color = "darkorange", alpha = 0.3, pch = 15) +
  geom_point(data = coord_muestras_g, size = 2,
             aes(x = NMDS1, y = NMDS2, 
                 color = generacion, shape = grupo_tratamiento)) +
  geom_text_repel(data = coord_muestras_g, size = 2,
                  aes(x = NMDS1, NMDS2, label = tratamiento))+
	 geom_text_repel(data = filter(coord_genero, genero_sub %in% abundancia_sub),
	                  aes(x = NMDS1, y = NMDS2, label = genero_sub), 
	                  color = "grey60", size = 2, min.segment.length = 100000) +
  scale_color_manual(values = paleta_gens) + 
  labs(subtitle = paste("Stress: ", round(nmds_genero$stress, 4), sep = "")) +
  scale_shape_manual(values = c("Control" = 17, "Tratamiento" = 19)) +
  theme_linedraw() 

ggsave("figuras/nmds_genero_1.png", dpi = 300, units = "cm", 
       height = 20, width = 20)


nmds_phylum <- ordinate(phylum_glom_obj, method = "NMDS",distance = "bray")

coord_phylum <- scores(nmds_phylum, display = "species") %>% data.frame() %>% 
	rownames_to_column(var = "OTU") %>% 
	left_join(tax, by = "OTU")

coord_muestras_p <- scores(nmds_phylum, display = "sites") %>% data.frame() %>% 
	rownames_to_column(var = "muestra") %>% 
	left_join(info, by = "muestra")


ggplot() + 
  #	geom_point(data = coord_phylum, aes(x = NMDS1, NMDS2), 
  #		   color = "darkorange", alpha = 0.5, pch = 15) +
  geom_text_repel(data = coord_phylum, aes(x = NMDS1, NMDS2, label = Phylum), 
                  size = 2, color = "grey60")+
  geom_point(data = coord_muestras_p, size = 2,
             aes(x = NMDS1, y = NMDS2, 
                 color = generacion, shape = grupo_tratamiento)) +
  geom_text_repel(data = coord_muestras_p, size = 2,
                  aes(x = NMDS1, NMDS2, label = tratamiento))+
  labs(subtitle = paste("Stress: ", round(nmds_phylum$stress, 4), sep = "")) +
  scale_color_manual(values = paleta_gens) + 
  scale_shape_manual(values = c("Control" = 17, "Tratamiento" = 19)) +
  theme_linedraw() 


ggsave("figuras/nmds_phylum.png", dpi = 300, units = "cm", 
       height = 20, width = 20)

```


## Enriquecimiento OTUs


### Controles vs tratamientos por generación

Agrupar las muestras por generación

```{r}

otus_tomate <- tomate_raw

sample_data(otus_tomate)$grupo_generacion <- as.character(sample_data(otus_tomate)$grupo_generacion)
otus_por_gen <- merge_samples(otus_tomate, group = "grupo_generacion")

```

Corregir metadatos del nuevo objeto

```{r}

sample_data(otus_por_gen)$grupo_generacion <- rownames(sample_data(otus_por_gen))
sample_data(otus_por_gen)

grupos_trat <- ifelse(str_detect(grupos, "C75"), "Control", "Tratamiento")
sample_data(otus_por_gen)$grupo_tratamiento <- grupos_trat

grupos <- sample_data(otus_por_gen)$grupo_generacion
sample_data(otus_por_gen)$tratamiento <- str_replace(grupos, ".+_", "")
sample_data(otus_por_gen)$generacion <- str_replace(grupos, "_.+", "")

sample_data(otus_por_gen)$replica <- NULL
sample_data(otus_por_gen)$grupo_generacion_muestra <- NULL

sample_data(otus_por_gen)

```

Correr DESeq de controles de cada generación contra los tratamientos, valores 
negativos corresponden a los OTUs sobrerepresentados en los controles

```{r}

deseq_controles_gen <- list()
for (gen in c("G1", "G2", "G3")){

  otus_gen <- subset_samples(otus_por_gen, generacion == gen)

  obj_deseq_gen <- phyloseq_to_deseq2(otus_gen, design = ~ grupo_tratamiento)
  obj_deseq_gen$grupo_tratamiento <- relevel(obj_deseq_gen$grupo_tratamiento, 
                                            ref = "Control")
  deseq_gen <- DESeq(obj_deseq_gen, test = "Wald", fitType = "local")
  out_gen <- results(deseq_gen) %>%
    as.data.frame() %>% rownames_to_column(var = "OTU") %>% tibble() %>% 
    mutate(generacion = gen)

  deseq_controles_gen[[gen]] <- out_gen

}

```

Gráficas de géneros sobrerepresentados en controles y tratamientos

```{r}

deseq_cont_df <- do.call(rbind, deseq_controles_gen)

#deseq_cont_df %>%
#  filter(padj < 0.05, log2FoldChange < 0) %>%
#  left_join(tax_sub, by = "OTU") %>%
#  select(OTU, tax_completa, genero_sub, baseMean:generacion) %>%
#  dplyr::count(genero_sub, generacion) %>%
#  arrange(generacion, desc(n)) %>%
#  write_tsv("tablas/generos_sob_controles.tsv")

tax_mod <- tax %>%
  mutate(categoria = ifelse(Phylum == "p__Pseudomonadota", Class, Phylum)) 

df_deseq_plot <- deseq_cont_df  %>% 
  filter(padj < 0.05) %>%
left_join(select(tax_sub, OTU, genero_sub), by ="OTU") %>%
left_join(tax_mod, by = "OTU") %>%
mutate(categoria = factor(categoria, levels = orden_categoria$categoria, 
                          ordered = TRUE))

df_deseq_plot %>%
  ggplot(aes(x = log2FoldChange,
             y = reorder(genero_sub, log2FoldChange),
             color = categoria)) +
  geom_point(size = 2) +
  geom_vline(xintercept = 0, color = "grey60") +
  labs(y = NULL, color = NULL) +
  facet_grid(generacion~"A", scales = "free_y", space = "free_y") +
  theme_linedraw() + 
  theme(axis.text.y = element_text(size = 10))

ggsave(filename = "figuras/sobrerepresentados_conts_trat.svg", 
       dpi = 300, units = "cm", height = k80, width = 22)

```

Tabla con OTUs sobrerepresentados por generación

```{r}

deseq_cont_df %>%
  filter(padj < 0.05) %>%
  left_join(tax_sub, by = "OTU") %>%
  select(OTU, tax_completa, genero_sub, generacion, baseMean:stat, padj) %>%
  write_tsv("tablas/deseq_controles_vs_tratamientos.tsv")  

deseq_cont_df %>%
  filter(padj < 0.05) %>%
  left_join(tax_sub, by = "OTU") %>%
  select(OTU, tax_completa, genero_sub, generacion, baseMean:stat, padj)  %>%
  group_by(generacion, genero_sub) %>%
  dplyr::count() %>%
  arrange(generacion, desc(n)) %>%
  write_tsv("tablas/deseq_controles_vs_tratamientos_generos.tsv")  

```



Lista de OTUs sobrerepresentados en controles y en tratamientos

```{r}

otus_sob_controles <- deseq_cont_df %>%
  filter(log2FoldChange < 0, padj < 0.05) %>% 
  pull(OTU) %>% unique()

tax_sub %>%
  filter(OTU %in% otus_sob_controles) %>%
  group_by(genero_sub) %>%
  dplyr::count() %>%
  arrange(desc(n)) %>%print(n= 10000)

otus_sob_trat <- deseq_cont_df %>%
  filter(log2FoldChange > 0, padj < 0.05) %>% 
  pull(OTU) %>% unique()

```

### Enriquecimiento secuencial por generación dentro de cada tratamiento

Las comparaciones siguen la siguiente tabla, comparando A vs B en cada 
tratamiento:

|A|B|
|:---:|:---:|
|S|G1|
|G1|G2|
|G2|G3|

Los valores de log2FoldChange negativos corresponden a los géneros
sobrerepresentados en A y los positivos a B.

```{r}

comparaciones <- c("G1" = "S", "G2" = "G1", "G3" = "G2")

deseq_multiple <- list()
for (k in c("S2", "S3", "S4")) {

  list_temp <- list()
  for (i in names(comparaciones)) {
    comp_temp <- paste(k,comparaciones[[i]], i, sep = "_")
    print(comp_temp)
    temp_sub <- subset_samples(tomate_raw, 
                               tratamiento == k & generacion %in% c(comparaciones[[i]], i))

    sample_data(temp_sub)$generacion <- factor(sample_data(temp_sub)$generacion, 
                                               ordered = FALSE)

    temp_deseq <- phyloseq_to_deseq2(temp_sub, design = ~ generacion)
    temp_deseq$generacion <- relevel(temp_deseq$generacion, ref = comparaciones[[i]])

    temp_analisis <- DESeq(temp_deseq, test = "Wald", fitType = "local")
    out_temp <- results(temp_analisis) %>%
      as.data.frame() %>% rownames_to_column(var = "OTU") %>% tibble()

    list_temp[[comp_temp]] <- out_temp
  }
  deseq_multiple[[k]] <- list_temp
}


lab_comparaciones <- c("SvsG1", "G1vsG2", "G2vsG3")

for (suelo in c("S2", "S3", "S4")) {
  df_temp <- data.frame()

  for (gen in 1:3) {
  df_suelo <- deseq_multiple[[suelo]][[gen]] %>%
    filter(padj < 0.05) %>%
    mutate(tratamiento = suelo, 
           comparacion = lab_comparaciones[gen]) %>%
    left_join(tax_sub, by = "OTU") %>%
    select(OTU, tax_completa, genero_sub, baseMean:comparacion)

  df_temp <- rbind(df_temp, df_suelo)
  }

  write_tsv(df_temp, 
            file = paste("tablas/deseq_por_generacion", suelo, ".tsv", sep = "" ))
}


```

## Microbioma núcleo

### Núcleo por  suelo

Loop para obtener el microbioma núcleo de cada tratamiento, considerando
ausencias en el inóculo incial (G0). Al usar el objeto con las muestras
agrupadas por tratamiento y generación, la condición para fomrar parte del
núcleo es presencia en al menos una de las muestras de cada generación. 

Antes de calcular el núcleo se quitan los OTUs que estuvieron representados en 
los controles con el vector `otus_sob_controles` .

```{r}

nucleos_tratamientos <- list()
for (suelo in c("S2", "S3", "S4")) {

  sub_suelo <- subset_samples(otus_por_gen, 
                              tratamiento == suelo)


  bin_sub_suelo <- ifelse(t(otu_table(sub_suelo)) > 0, 1, 0) %>% 
    data.frame() %>%
    filter(!rownames(.) %in% otus_sob_controles)

  # Otus ausentes en los suelos, pero presentes en las raíces

  sin_suelo <- bin_sub_suelo[rowSums(bin_sub_suelo) == 3,] %>% data.frame() %>%
                    filter(if_all(starts_with(c("G1", "G2", "G3")), ~ .x == 1)) %>%
                    rownames()

  con_suelo <- bin_sub_suelo[rowSums(bin_sub_suelo) == 4,] %>%
    data.frame() %>% rownames()

  if (sum(sin_suelo %in% con_suelo) == 0 & sum(con_suelo %in% sin_suelo) == 0) {

    suelo_nucleo <- c(sin_suelo, con_suelo)
    nucleos_tratamientos[[suelo]] <- suelo_nucleo

    png(filename = paste("figuras/upset_", suelo,".png", 
                         sep = ""), units = "cm", 
        width = 12, height = 8, res = 300)
    plot_suelo <- upset(data = bin_sub_suelo,
                        sets = colnames(bin_sub_suelo), nintersects = 20,
                        point.size = 2, line.size = 1, 
                        mainbar.y.label = "OTUs compartidos", 
                        sets.x.label = "OTUs por condicion",
                        decreasing = TRUE, order.by = "degree")
    print(plot_suelo)
    dev.off()
  }

}

nucleos_tratamientos

lapply(nucleos_tratamientos, length)

```

###  Núcleo de las tres generaciones

Crear nuevo objeto con los núcleos de cada tratamiento para calcualr la 
intersección entre todos los tratamientos de las tres generaciones

```{r}

otus_tratamientos <- unlist(nucleos_tratamientos) %>% unname() %>%unique()

bin_tratamientos <- data.frame(otus_tratamientos,
       "S2" = as.numeric(otus_tratamientos %in% nucleos_tratamientos[[1]]),
       "S3" = as.numeric(otus_tratamientos %in% nucleos_tratamientos[[2]]),
       "S4" = as.numeric(otus_tratamientos %in% nucleos_tratamientos[[3]])) %>%
column_to_rownames(var = "otus_tratamientos")


png(filename = "figuras/upset_tratamientos.png", 
    units = "cm", width = 20, height = 10, res = 300)
upset(data = bin_tratamientos,
      sets = colnames(bin_tratamientos), nintersects = 20,
      point.size = 2, line.size = 1, 
      mainbar.y.label = "OTUs compartidos", 
      sets.x.label = "OTUs por condicion",
      decreasing = TRUE, order.by = "degree")
dev.off()

otus_nucleo_tratamientos <- bin_tratamientos %>%
  filter(rowSums(.) == 3) %>% rownames()

tax_sub %>%
  filter(OTU %in% otus_nucleo_tratamientos) %>%
  group_by(genero_sub) %>%
  dplyr::count() %>%
  arrange(desc(n)) %>%
  ungroup() %>%
  write_tsv("tablas/generos_nucleo.tsv")


tax_sub %>%
  filter(OTU %in% otus_nucleo_tratamientos)

  write_tsv("tablas/OTUs_nucleo.tsv")


nucleo_tratamientos <- bin_tratamientos %>% 
  mutate(rsum = rowSums(.)) %>% 
  filter(rsum == 3) %>% 
  rownames()

```


## Gráficas abundancia y enriquecimiento por suelo

Graficar mapa de calor de los géneros del nucleo de cada suelo, junto con
información de DESeq generacional y sobrerepresentado entre controles

### Datos para cada grafica

```{r}


dfs_graficas <- list()
for (suelo in c("S2", "S3", "S4")) {

  orden_suelo <- tabla_ab_rel_otus  %>% 
    filter(tratamiento == suelo) %>% 
    group_by(OTU) %>% 
    summarise(prom_ab = mean(ab_rel)) %>% 
    arrange(desc(prom_ab))

  # Generar data frame con  OTUs que tuvieron abundancia diferencial en alguna 
  # de las generaciones

  dif_suelo <- data.frame()
    for (i in names(comparaciones)) {
      comp_temp <- paste(suelo,comparaciones[[i]], i, sep = "_")

      dif_suelo <- rbind(dif_suelo, 
        deseq_multiple[[suelo]][[comp_temp]] %>%
          filter(padj <= 0.05) %>%
          mutate(comp = comp_temp))
    }

  dif_suelo <- dif_suelo  %>%
    left_join(tax_sub, by = "OTU") %>%
    mutate(comp = factor(comp, 
                         levels = paste(suelo, c("_S_G1", "_G1_G2", "_G2_G3"), sep = ""), 
                         ordered = TRUE))


  # Crear vectores para dividir OTUs en persistentes(nucleo) y transitorios 
  suelo_nucleo <- nucleos_tratamientos[[suelo]] 

  deseq_extendido <- dif_suelo %>% pull(OTU) %>%setdiff(suelo_nucleo) %>%
    setdiff(otus_sob_controles)

  # Generar vector con OTUs enriquecidos para cada tratamiento

  df_base <- tabla_ab_rel_otus %>% 
    filter(tratamiento == suelo) %>%
    select(OTU, ab_rel, muestra, generacion, grupo_replica, tax_completa, genero_sub) %>%
    left_join(orden_suelo, by = "OTU") %>%
    left_join(select(dif_suelo, OTU, baseMean:stat, comp, padj), by = "OTU") %>%
    filter(OTU %in% c(deseq_extendido, suelo_nucleo)) %>%
    mutate(bin = ifelse(OTU %in% otus_sob_trat, "1", "0"), 
           tipo = factor(ifelse(OTU %in% suelo_nucleo, "core", "ext"), 
                         levels = c("core", "ext"), ordered  = TRUE), 
           OTU_sub = paste(OTU, genero_sub, sep = "_")) %>%
    replace_na(list(log2FoldChange = 0))

    dfs_graficas[[suelo]] <- df_base

}


for (gen in c("S2", "S3", "S4")){

  dfs_graficas[[gen]] %>%
    separate_wider_delim(cols = comp, delim = "_", 
                         names = c("trat", "ref", "sob_en")) %>% 
    select(OTU, tax_completa, genero_sub, tipo, ab_rel, generacion,
         grupo_replica, sob_en, baseMean:stat, padj) %>%
    arrange(generacion, desc(padj)) %>%
    mutate(tipo = ifelse(tipo == "core", "persistente", "transitorio")) %>%
    write_tsv(paste("tablas/persistentes_y_transitorios_", gen, ".tsv", sep = ""))

}

```

### Gráficas de enriquecimineto por generación, separadas por tratamiento

```{r}

for (suelo in c("S2", "S3", "S4")) {

  core_plot <- dfs_graficas[[suelo]] %>%
    ggplot(aes(x = grupo_replica, y = reorder(OTU_sub, prom_ab), 
               fill = ab_rel)) +
    geom_raster() + 
    scale_fill_gradient(low = "#e4e4e4ff", high = "black", 
                        na.value = "steelblue4", trans = "log10", 
                        guide = guide_colorbar(barwidth = unit(.4, "cm"),
                                               barheight = unit(2, "cm"))) +
    labs(y = NULL, x = NULL, fill = "Abundancia\nrelativa", 
         title = suelo) + 
    scale_y_discrete(labels = NULL) +
    scale_x_discrete(labels = NULL) +
    facet_grid(tipo~generacion, scales = "free", space = "free") +
    theme_linedraw() +
    theme(strip.text.y = element_blank(), 
          strip.background.y = element_blank(), 
          axis.text.x = element_text(angle = -270),
          axis.text = element_text(size = 5), 
          axis.ticks = element_blank(), 
          legend.key.size = unit(1, "cm"), 
          legend.text = element_text(size = 5), 
          legend.title = element_text(size = 5),
          panel.grid = element_blank())


  deseq_plot <- dfs_graficas[[suelo]] %>%
    ggplot(aes(x = comp, y = reorder(OTU, prom_ab), 
               fill = log2FoldChange)) +
    geom_raster(na.rm = TRUE) +
    labs(x = NULL, y = NULL) +
    scale_y_discrete(labels = NULL) +
    scale_x_discrete(labels = NULL)+
    scale_fill_gradient2(low = "#600076ff", high = "#338b00ff", 
                         mid = "white",
                         guide = guide_colorbar(barwidth = unit(.4, "cm"),
                                                barheight = unit(2, "cm"))) +
    theme_linedraw() +
    facet_grid(tipo~"A", scales = "free", space = "free") +
    theme(strip.text = element_blank(), 
          strip.background = element_blank(), 
          axis.text.x = element_text(angle = -270),
          axis.text = element_text(size = 5), 
          axis.ticks = element_blank(), 
          legend.key.size = unit(1, "cm"), 
          legend.text = element_text(size = 5), 
          legend.title = element_text(size = 5),
          panel.grid = element_blank())

  core_plot + deseq_plot + 
    plot_layout(nrow = 1, guides = "collect", 
                widths = c(7, 3))

  ggsave(filename = paste("figuras/otus_core_ext_", suelo, "_compacto.png", sep = ""), 
         dpi = 300, units = "cm", height = 15, width = 12, limitsize = FALSE)

  ggsave(filename = paste("figuras/otus_core_ext_", suelo, "_compacto.svg", sep = ""), 
         dpi = 300, units = "cm", height = 15, width = 12, limitsize = FALSE)

}

```

## Abundancia relativa de OTUS sobrerepresentados en la comparación entre generaciones

Data frame con los OTUs enriquecidos por generación y tratamiento

```{r}

div_otus_sob <- data.frame()
for (suelo in c("S2", "S3", "S4")) {
  lista_gtratOen <- c()
  for (gen in 1:3) {
    otus_suelo <-  deseq_multiple[[suelo]][[gen]] %>%
      filter(padj < 0.05, log2FoldChange > 0) %>% 
      mutate(sob_gen = paste("G", gen, sep = ""), 
             sob_trat = suelo)
      div_otus_sob <- rbind(otus_suelo, div_otus_sob)
  }
}

```

Tabla con los datos completos de OTUs enriquecidos

```{r}

sub_otus_sob <- div_otus_sob %>% pull(OTU) %>% unique()

tabla_ab_rel_otus %>%
  filter(OTU %in% sub_otus_sob, 
         tratamiento != "C75") %>%
left_join(orden_ab_otus, by = "OTU") %>%
left_join(select(div_otus_sob, OTU, sob_gen, sob_trat,baseMean:stat, padj), by = "OTU") %>%
select(OTU, tax_completa, genero_sub, generacion, ab_rel,
       tratamiento, replica, sob_trat, sob_gen, 
       baseMean:stat, padj) %>%
arrange(sob_gen, desc(padj)) %>%
write_tsv("tablas/otus_enriquecidos_generaciones.tsv")

```

Gráfica ordenada por abundancia relativa y separando los OTUs por 
la generación en donde ocurrió el enriquecimiento

```{r}

tabla_ab_rel_otus %>%
  filter(OTU %in% sub_otus_sob, 
         tratamiento != "C75") %>% 
  left_join(orden_ab_otus, by = "OTU") %>%
  left_join(div_otus_sob, by = "OTU") %>% 
  mutate(grupo_generacion_muestra = factor(grupo_generacion_muestra, 
                                           levels = orden_grupos, 
                                           ordered = TRUE)) %>% pull(grupo_generacion_muestra)

  ggplot(aes(x = grupo_generacion_muestra,
             y = reorder(OTU, prom_ab),
             fill = ab_rel)) + 
  geom_raster() + 
  scale_fill_gradient(low = "#e4e4e4ff", high = "#052629ff", 
                      na.value = "steelblue4", trans = "log10") +
  facet_grid(sob_gen~tratamiento, scales = "free", space = "free") +
  scale_y_discrete(labels = NULL)+
  #  scale_x_discrete(labels = NULL)+
  labs(y = NULL, x = NULL, fill = "Abundancia\nrelativa") +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = -270, size = 6),
        axis.text = element_text(size = 5), 
        axis.ticks = element_blank(), 
        legend.key.size = unit(1, "cm"), 
        panel.grid = element_blank())

ggsave(filename = "figuras/otus_enriquecidos_ordenados.png", 
       width = 12, height = 15, units = "cm", dpi = 300)

ggsave(filename = "figuras/otus_enriquecidos_ordenados.svg", 
       width = 12, height = 15, units = "cm", dpi = 300)

```

Conteo de OTUs enriquecidos por tratamientos, por generación

```{r}

conteos_sob_gen <- data.frame("V1" = 1:3)

for (gen in 1:3) {
  lista_gen <- c()
  for (suelo in c("S2", "S3", "S4")) {
    otus_suelo <-  deseq_multiple[[suelo]][[gen]] %>%
      filter(padj < 0.05, log2FoldChange > 0) %>% pull(OTU) %>%
      length()
    lista_gen <- c(lista_gen, otus_suelo)
  }
  conteos_sob_gen[gen] <- lista_gen
}

colnames(conteos_sob_gen) <- c("G1", "G2", "G3")
conteos_sob_gen$tratamiento <- c(c("S2", "S3", "S4"))

conteos_sob_gen %>%
  select(tratamiento, G1:G3) %>%
  write_tsv("tablas/conteos_sob_generaciones.tsv")

```

### Intersección de OTUs del núcleo de las tres generaciones y OTUs sobrerepresentados en tratamientos respecto al control

Tabla de OTUs y géneros del núcleo sobrerepresentados

```{r}

tax_sub %>%
  filter(OTU %in% sob_nucleo) %>% 
  left_join(filter(deseq_cont_df, padj < 0.05), 
            by = "OTU")  %>%
  select(OTU, tax_completa, genero_sub, generacion, padj, baseMean:stat, padj) %>%
  arrange(generacion, desc(padj)) %>%
  write_tsv("tablas/otus_nucleo_sobrerepresentados.tsv")

## Tabla con géneros
tax_sub %>%
  filter(OTU %in% sob_nucleo) %>% 
  left_join(filter(deseq_cont_df, padj < 0.05), 
            by = "OTU")  %>%
  select(genero_sub, OTU,generacion, padj,) %>%
  pivot_wider(names_from = generacion, values_from = padj) %>%
  arrange(genero_sub) %>%
  write_tsv("tablas/generos_nucleo_sob_signif.tsv")

```

Gráfica de abundancia relativa y mapa de calor mostrando log2FoldChange 

```{r}

sob_nucleo <- otus_nucleo_tratamientos[otus_nucleo_tratamientos %in% otus_sob_trat]

orden_sob_nucleo <- tabla_ab_rel_otus %>%
  filter(OTU %in% sob_nucleo) %>%
  group_by(OTU) %>%
  summarise(prom_ab = mean(ab_rel)) %>%
  arrange(desc(prom_ab))

ab_rel_nucleo_sob <- tabla_ab_rel_otus %>%
  filter(OTU %in% sob_nucleo, 
         tratamiento != "C75") %>%
  left_join(orden_sob_nucleo, by = "OTU") %>%
  mutate(OTU_sub = paste(OTU, genero_sub, sep = "_"), 
         grupo_generacion_muestra = factor(grupo_generacion_muestra, 
                                           orden_grupos, 
                                           ordered  = TRUE)) %>%
  ggplot(aes(x = grupo_generacion_muestra, 
             y = reorder(OTU_sub, prom_ab), 
             fill = ab_rel)) +
  geom_raster()  +
  labs(x = NULL, y = NULL, fill = "Abundancia\nrelativa") +
  scale_fill_gradient(low = "#e4e4e4ff", high = "#052629ff", 
                      na.value = "steelblue4", trans = "log10") +
  facet_grid(~tratamiento, scales = "free", space = "free") + 
  theme_linedraw() + 
  theme(axis.text.x = element_text(angle = -270), 
        axis.text = element_text(size = 6))


deseq_nucleo_sob <- tax_sub %>%
  filter(OTU %in% sob_nucleo) %>% 
  left_join(filter(deseq_cont_df, padj < 0.05),
            by = "OTU") %>%
  left_join(orden_sob_nucleo, by = "OTU") %>%
  ggplot(aes(x = generacion, y = reorder(OTU, prom_ab), 
             fill = log2FoldChange)) + 
  geom_raster() +
  labs(x = NULL, y = NULL, fill = "log2FoldChange") + 
	scale_fill_gradient2(low = "purple2",mid = "white", high = "palegreen3") + 
  theme_linedraw() + 
  theme(strip.text = element_blank(), 
        strip.background = element_blank(), 
        axis.text.y = element_blank(), 
        axis.text.x = element_text(angle = -270), 
        axis.ticks = element_blank(), 
        legend.key.size = unit(1, "cm"), 
        panel.grid = element_blank())

ab_rel_nucleo_sob + deseq_nucleo_sob +
    plot_layout(nrow = 1, guides = "collect", 
                widths = c(7, 1))

ggsave(filename = "figuras/heatmap_core_tratamientos_sob.png", 
       units = "cm", height = 20, width = 20, dpi = 300)

ggsave(filename = "figuras/heatmap_core_tratamientos_sob.svg", 
       units = "cm", height = 20, width = 20, dpi = 300)

```




## Fenotipo vegetal y pH


### PCA solo de variables fenotipo y pH

```{r}

tabla_pca <- fenotipo %>%
  group_by(generacion, tratamiento, replica) %>%
  summarise(prom_diametro = mean(diametro_tallo_rel, na.rm = TRUE), 
            prom_longitud = mean(longitud_tallo_rel, na.rm = TRUE), 
            prom_peso_aereo = mean(peso_seco_aereo_rel, na.rm = TRUE), 
            prom_peso_raiz = mean(peso_seco_raiz_rel, na.rm = TRUE),
            prom_pH = mean(pH, na.rm = TRUE)) %>%
  ungroup()


pca_fenotipo <-  prcomp(x = tabla_pca[,-(1:3)], scale = TRUE)

scores <- cbind(tabla_pca[,1:3], data.frame(pca_fenotipo$x)) %>%
  mutate(grupo_tratamiento = ifelse(tratamiento == "C75", "Control", "Tratamiento"), 
         tratamiento_label = ifelse(tratamiento == "C75", "", tratamiento))

loadings <- rownames_to_column(data.frame(pca_fenotipo$rotation), var = "variable") %>%
  mutate(var_lab = c("Diámetro tallo", "Longitud tallo", 
                     "Peso seco aéreo", 
                     "Peso seco raíz"))


var_componentes <- round((pca_fenotipo$sdev**2 / sum(pca_fenotipo$sdev**2) * 100), 2)

ggplot() + 
  geom_point(data = scores, aes(x = PC1, y = PC2, 
                                shape = grupo_tratamiento, 
                                color = generacion), 
             size = 3) + 
  geom_segment(data = loadings, aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_text(data = loadings, aes(x = PC1, y = PC2, label = var_lab),
                  size = 2) + 
  geom_text_repel(data = scores, aes(x = PC1, y = PC2, label = tratamiento_label), 
                  size = 2, color = "grey60") +
  scale_color_manual(values = paleta_gens) + 
  scale_shape_manual(values = c(17, 19)) + 
  labs(x = paste("PC1 [", var_componentes[1], "%]", sep = ""), 
       y = paste("PC2 [", var_componentes[2], "%]", sep = ""), 
       color = "Generación", shape = NULL) +
  theme_linedraw()


ggsave("figuras/pca_fenotipo.png", units = "cm", 
       height = 12, width = 16, dpi = 300)
ggsave("figuras/pca_fenotipo.svg", units = "cm", 
       height = 12, width = 16, dpi = 300)


```

### Pruebas estadisticas fenotipo 


```{r}

vars_fenotipo <- colnames(fenotipo)[16:20]

anova_fenotipo <- list()
for (var in vars_fenotipo){

  lista_var <- list()
  for (gen in c("G1", "G2", "G3")) {
    lista_gen <- list()

    data_gen <- fenotipo %>%
      filter(generacion == gen)

    anova_var <- aov(formula = get(var) ~ tratamiento, 
                     data = data_gen)

    normalidad <- shapiro.test(anova_var$residuals)
    varianza <- bartlett.test(formula = get(var) ~ tratamiento, 
                              data = data_gen)

    lista_gen[[1]] <- normalidad
    lista_gen[[2]] <- varianza
    lista_gen[[3]] <- summary(anova_var)
    post_hoc <- TukeyHSD(anova_var)
    lista_gen[[4]] <- post_hoc

    lista_var[[gen]] <- lista_gen

  }

  anova_fenotipo[[var]] <- lista_var
}

anova_fenotipo


anova_fenotipo[[2]][[4]]$grupo_generacion %>%
  data.frame() %>% 
  filter(p.adj < 0.05)


kruskal_fenotipo <- list()
for (var in vars_fenotipo){
  print(var)
  lista_var <- list()
  kruskal_var <- kruskal.test(formula = get(var) ~ grupo_generacion, 
                                 data = fenotipo)

  lista_var[[1]] <- kruskal_var

  if (kruskal_var$p.value < 0.05) {
    wilcox_par <- pairwise.wilcox.test(x = unlist(fenotipo[,var]), 
                         g = fenotipo$grupo_generacion)
    lista_var[[2]] <- wilcox_par
  }

    kruskal_fenotipo[[var]] <- lista_var
}

```

### Gráficas de fenotipo


```{r}

long_tallo <- fenotipo %>%
  ggplot(aes(x = generacion, y = longitud_tallo_rel, color = tratamiento)) + 
  geom_boxplot() + 
  geom_point(position = position_dodge(width = 0.75)) + 
  geom_signif(y_position = c(1, .9, .8), 
              xmin = c(.72, 1.14, 1.72), 
              xmax = c(1.28, 1.28, 2.28), 
              annotation = c("**", "*", "*"), color = "black") +
  scale_color_manual(values = paleta, labels = dict_tratamientos) +
  labs(x = NULL, y = "Largo de tallo relativo (cm)",
       color = NULL, title = "A") +
  theme_linedraw()

diam_tallo <- fenotipo %>%
  ggplot(aes(x = generacion, y = diametro_tallo_rel, color = tratamiento)) + 
  geom_boxplot() + 
  geom_point(position = position_dodge(width = 0.75)) + 
  scale_color_manual(values = paleta, labels = dict_tratamientos) +
  labs(x = NULL, y = "Diámetro de tallo relativo (mm)",
       color = NULL, title = "B") +
  theme_linedraw()

long_tallo + diam_tallo + plot_layout(guides = "collect", ncol = 1)

ggsave("figuras/fenotipo_1.png", dpi = 300, units = "cm",
       height = 16, width = 16)

## Peso seco 

peso_total <- fenotipo %>%
  ggplot(aes(x = generacion, y = peso_seco_total_rel, 
             color = tratamiento)) + 
  geom_boxplot() + 
  geom_point(position = position_dodge(width = 0.75)) + 
  scale_color_manual(values = paleta, labels = dict_tratamientos) +
  labs(x = NULL, y = "Peso seco (g)",
       color = NULL, title = "A", subtitle = "Total") +
  theme_linedraw()

peso_tallo <- fenotipo %>%
  ggplot(aes(x = generacion, y = peso_seco_aereo_rel,
             color = tratamiento)) + 
  geom_boxplot() + 
  geom_point(position = position_dodge(width = 0.75)) + 
  scale_color_manual(values = paleta, labels = dict_tratamientos) +
  labs(x = NULL, y = "Peso seco (g)",
       color = NULL, title = "B", subtitle = "Aéreo") +
  theme_linedraw()

peso_raiz <- fenotipo %>%
  ggplot(aes(x = generacion, y = peso_seco_raiz_rel, color = tratamiento)) + 
  geom_boxplot() + 
  geom_point(position = position_dodge(width = 0.75)) + 
  scale_color_manual(values = paleta, labels = dict_tratamientos) +
  labs(x = NULL, y = "Peso seco (g)",
       color = NULL, title = "C", subtitle = "Raíz") +
  theme_linedraw()

peso_total + peso_tallo + peso_raiz + 
  plot_layout(guides = "collect", ncol = 1)

ggsave("figuras/fenotipo_2.png", dpi = 300, units = "cm",
       height = 18, width = 16)

fenotipo %>%
  ggplot(aes(x = peso_seco_aereo_rel, y = peso_seco_raiz_rel, 
             color = tratamiento)) + 
  geom_point() +
  facet_wrap(~generacion, nrow = 1) + 
  scale_color_manual(values = paleta)


ggsave("figuras/corr_peso_seco.png", dpi = 300, units = "cm",
       height = 10, width = 20)
```

### Gráfica de pH

```{r}

fenotipo %>%
  ggplot(aes(x = tratamiento, y = pH, color = generacion)) +
  geom_boxplot() +
  geom_point(position = position_dodge(width = 0.75)) +
  labs(color = "Generación", x = NULL) + 
  scale_color_manual(values = paleta_gens) +
  theme_linedraw()

ggsave("figuras/pH.png", units = "cm", 
       width = 12, height = 10, dpi = 300)

```

## Infecciones con B. cinerea

```{r}

hojas <- read_csv("infecciones.csv", col_types = "cccfffd") %>%
  filter(tratamiento != "C100") %>%
  mutate(tratamiento = factor(tratamiento,
                              levels = c("C75", "S2", "S3", "S4"),
                              ordered = TRUE))

```

### Area de infeccion
Los datos van a estar agrupados por replica, promediando 
los valores de las cuatro hojas


```{r}

hojas_G2 <- hojas %>% filter(generacion == "G2")
hojas_G3 <- hojas %>% filter(generacion == "G3")

modelo_infeccion_G2 <- aov(formula = area_infeccion ~ tratamiento, data = hojas_G2)
shapiro.test(modelo_infeccion_G2$residuals)
bartlett.test(formula = area_infeccion ~ tratamiento, data = hojas_G2)

modelo_infeccion_G3 <- aov(formula = area_infeccion ~ tratamiento, data = hojas_G3)
shapiro.test(modelo_infeccion_G3$residuals)
bartlett.test(formula = area_infeccion ~ tratamiento, data = hojas_G3)

TukeyHSD(modelo_infeccion_G2) %>%
   tidy() %>%
   filter(adj.p.value <= 0.05)

# No hay diferencias significativas entre tratamientos en G3
kruskal.test(formula = area_infeccion ~ tratamiento, data = hojas_G3)


hojas %>%
  ggplot(aes(x = generacion, y = area_infeccion, color = tratamiento)) +
  geom_boxplot() +
  geom_point(position = position_dodge(width = .75)) +
  geom_signif(y_position = c(2, 1.75, 1.5),
              xmin = c(.75, .9, 1.1), xmax = c(1.28, 1.28, 1.28),
              annotation = c("***", "***", "*"),
              color = "black", tip_length = .02) +
  scale_color_manual(values = paleta, labels = dict_tratamientos) +
  scale_x_discrete(labels = c("G2", "G3")) +
  labs(x = NULL, y = "Área de infección (cm^2)", color = NULL) +
  theme_bw()

ggsave("figuras/infecciones.png", dpi = 300, units = "cm",
       height = 8, width = 18)

```

# Mapa de sitios de colecta

```{r}

limite_mexico <- read_sf("limite_mexico/contdv250kgw.shp")

localidades <- data.frame(localidad = c("san_cayetano", 
                                        "Tziritzicuaro","Guanajuato"),
                          clave = c("S2", "S3", "S4"),
                          lat = c(19.38990639,
                                  19.95783333,
                                  19.78004806
                          ), 
                          long = c(-99.70802806,
                                   -100.5026944,
                                   -100.677875
                          ))

ggplot(data = localidades, aes(x = long, y = lat, color = clave))+
  geom_sf(data = limite_mexico, aes(geometry = geometry), 
          fill = "gray90", color = "gray20", inherit.aes = FALSE, 
          linewidth = 1/3) +
  geom_point(show.legend = FALSE) +
  geom_text(aes(label = clave), 
            nudge_x = c(0.7, 0.5, -0.5),
            nudge_y = c(0, 0.5, -0.5), show.legend = FALSE ) +
  labs(x = NULL, y = NULL, color = NULL) +
  scale_color_manual(values = paleta) +
  theme_linedraw() +
  theme(panel.background = element_rect(fill = "#ecf2fbff"))

ggsave(filename = "mapa_colecta_suelos.png", units = "cm", 
       height = 10, width = 15, dpi = 300)

```

