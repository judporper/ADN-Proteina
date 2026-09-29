# Del ADN a la proteína

**Extensiones con Biopython**

Universidad de Las Palmas de Gran Canaria · Grado en Ciencia e Ingeniería de Datos · Bioinformática
Sofía Travieso García, Judith Portero Pérez · 29 de septiembre de 2026

## 1. Introducción y Objetivos

Este trabajo recorre el dogma central de la biología molecular (replicación, transcripción y traducción), el splicing alternativo y la relación entre secuencia y estructura de las proteínas. Cada ejercicio tiene dos partes: una **parte manual**, razonada a mano, y una **extensión con Biopython** que la reproduce y amplía en un notebook (`Del_ADN_a_la_Proteína.ipynb`) con datos reales: el gen *lacZ* de *E. coli*, las isoformas de FGFR2, la hemoglobina humana (PDB 1GZX) y el gen TP53.

Repositorio: https://github.com/judporper/ADN-Proteina

## 2. Métodos y Datos

Se usó Python 3 y Biopython en un notebook de Jupyter. Los datos son públicos y están en el repositorio en formato FASTA.

| Ej. | Datos | Origen |
|---|---|---|
| 2 | *lacZ* de *E. coli* K-12 (NC_000913.3:c366305-363231, 3075 nt) | NCBI |
| 4 | CDS de los 46 transcritos de FGFR2 (ENSG00000066468) | Ensembl |
| 5 | Hemoglobina humana, cadenas α y β (1GZX) | RCSB PDB |
| 6 | TP53, transcrito canónico ENST00000269305.9 (1182 nt) | Ensembl |

| Función | Qué hace y por qué se usa |
|---|---|
| `complement()`, `reverse_complement()` | Empareja A–T y C–G; la segunda escribe la hebra en 5'→3', como se anota una secuencia. |
| `transcribe()` | Cambia T por U. Solo es correcta si recibe la hebra codificante en 5'→3'. |
| `complement_rna()` | Complementa la hebra molde con bases de ARN, como la ARN polimerasa. |
| `translate(to_stop=True)` | Aplica el código genético estándar y se detiene en el primer codón de paro. |
| `PairwiseAligner` | Alinea proteínas (BLOSUM62) para comparar isoformas. |
| `ProteinAnalysis`, `Bio.PDB`, `ShrakeRupley` | Composición de la secuencia, estructura secundaria y superficie accesible. |

## 3. Desarrollo de los Ejercicios

### 3.1. Ejercicio 1. Replicación del ADN

**Parte manual.** La duplicación es **semiconservativa**: cada hebra parental sirve de molde y las dos moléculas hijas conservan una hebra vieja y contienen una nueva. Para el fragmento `5'-ATG CCG TTA GCT-3'` / `3'-TAC GGC AAT CGA-5'` se obtiene:

```
Hija 1:  5'-ATG CCG TTA GCT-3'  (parental)      Hija 2:  5'-ATG CCG TTA GCT-3'  (nueva)
         3'-TAC GGC AAT CGA-5'  (nueva)                  3'-TAC GGC AAT CGA-5'  (parental)
```

- **Helicasa:** abre la doble hélice rompiendo los puentes de hidrógeno entre bases y forma la horquilla de replicación.
- **Primasa:** sintetiza un cebador corto de ARN que aporta el extremo 3'-OH que la ADN polimerasa necesita, porque esta no puede iniciar una cadena desde cero.
- **ADN polimerasa:** añade desoxirribonucleótidos en sentido 5'→3' leyendo el molde en 3'→5' (hebra líder de forma continua, hebra retardada en fragmentos de Okazaki) y corrige errores con su actividad exonucleasa 3'→5'.
- **Ligasa:** sella con enlaces fosfodiéster los huecos entre fragmentos de Okazaki, tras retirar los cebadores.

**Reflexión.** Si la polimerasa se equivoca y el error no se corrige, la base incorrecta queda en la hebra nueva y en la siguiente replicación se copia como si fuera correcta: es una mutación heredable. Según dónde caiga, puede ser silenciosa, cambiar un aminoácido o crear un codón de paro.

**Extensión con Biopython.** `complement()` devolvió `TACGGCAATCGA`, idéntica a la hebra obtenida a mano, y `reverse_complement()` la misma hebra escrita en 5'→3' (`AGCTAACGGCAT`).

### 3.2. Ejercicio 2. Transcripción del ADN a ARN

**Parte manual.** En `5'-ATG CCT GAA TGC-3'` / `3'-TAC GGA CTT ACG-5'`, la **cadena molde** es la inferior (3'→5'): la ARN polimerasa la lee en ese sentido y sintetiza el ARN en 5'→3' con U en lugar de T.

```
ADN codificante  5'-ATG CCT GAA TGC-3'
ADN molde        3'-TAC GGA CTT ACG-5'
ARNm             5'-AUG CCU GAA UGC-3'
```

El ARNm coincide con la hebra codificante con T→U. La **región promotora** no está en el fragmento: se encuentra en 5' del inicio de la transcripción, donde se une la ARN polimerasa, y no se transcribe. La **región codificante** es el propio fragmento, que empieza en el codón de inicio ATG (Met-Pro-Glu-Cys).

**Extensión con Biopython.** Se leyó *lacZ* desde FASTA (3075 nt) y se transcribió: empieza por AUG y se traduce a 1024 aa, la longitud de la β-galactosidasa real. Al cambiar la orientación (experimento pedido), la reversa complementaria da un ARNm que empieza por UUA y la complementaria sin invertir por UAC; en ambos casos la proteína ya no es la β-galactosidasa, porque `transcribe()` no orienta la hebra y eso lo decide quien se la pasa. Transcribir desde la hebra molde con `complement_rna()` da el mismo ARNm que `transcribe()` sobre la codificante.

### 3.3. Ejercicio 3. Traducción del ARNm a proteína

**Parte manual.** Para `5'-AUG UAU GCU UAA-3'`, el codón de inicio es AUG y el de paro UAA:

| Codón | AUG | UAU | GCU | UAA |
|---|---|---|---|---|
| Significado | inicio (Met) | Tyr | Ala | paro |

La cadena resultante es **Met-Tyr-Ala**.

**Reflexión.** Si el codón de inicio muta de AUG a GUG (Val), el ribosoma eucariota no reconoce el inicio y no traduce ese marco, salvo que haya otro AUG más adelante; en procariotas GUG puede iniciar la traducción, aunque con menor eficiencia. Si el codón de paro desaparece, la traducción continúa (*readthrough*) hasta encontrar otro paro en la secuencia siguiente y da una proteína más larga, con aminoácidos añadidos, probablemente inestable o no funcional.

**Extensión con Biopython.** `translate(to_stop=True)` da `MYA`, igual que a mano; `translate()` sin más lo marca como `MYA*`.

### 3.4. Ejercicio 4. Splicing alternativo

**Parte manual.** Con un gen de cinco exones se proponen tres combinaciones de ARNm maduro:

| Isoforma | Exones | Diferencia esperada en la proteína |
|---|---|---|
| A (completa) | 1-2-3-4-5 | Proteína con todos los dominios. |
| B (salto de un exón) | 1-2-4-5 | Pierde el segmento del exón 3. Si su longitud no es múltiplo de 3, cambia el marco de lectura y aparece un paro prematuro. |
| C (salto de dos exones) | 1-3-5 | Proteína más corta, sin los segmentos de los exones 2 y 4. |

**Reflexión.** Un solo gen produce varios ARNm maduros según los exones que se unan, y cada uno codifica una proteína distinta; con los exones 1 y 5 fijos, los tres intermedios ya permiten hasta 2³ = 8 combinaciones. Así aumenta la diversidad del proteoma sin aumentar el número de genes.

**Extensión con Biopython / Ensembl.** Se descargaron los CDS de los 46 transcritos de FGFR2 y se tradujeron: hay 39 longitudes de CDS distintas y proteínas de 51 a 822 aa. El notebook marca los CDS que no son múltiplo de 3 (anotaciones probablemente parciales), cuenta las proteínas realmente distintas (varios transcritos pueden dar la misma) y alinea cada una con la más larga. Se esperan dos tipos de diferencia: las isoformas más cortas pierden segmentos completos (dominios extracelulares, región transmembrana o dominio quinasa), y sin dominio quinasa o sin región transmembrana no pueden transmitir la señal del receptor; las de igual longitud con sustituciones concentradas corresponden a exones mutuamente excluyentes, como FGFR2-IIIb y IIIc, descritas en la literatura, que reconocen ligandos FGF distintos.

### 3.5. Ejercicio 5. Introducción a las proteínas

**Parte manual.** En `Met-Ile-Ser-Gly-Val-Lys-His`, el extremo **N** es Met (grupo amino libre) y el extremo **C** es His (grupo carboxilo libre), porque la cadena se sintetiza desde el extremo N hacia el C.

**Reflexión.** El orden de los aminoácidos (estructura primaria) determina cómo se pliega la proteína: los residuos hidrofóbicos tienden a agruparse en el interior, lejos del agua, y los polares en la superficie. Si una mutación cambia un residuo hidrofóbico interno por uno hidrofílico, se introduce polaridad en el núcleo, se rompe el empaquetamiento y el plegamiento puede desestabilizarse o perderse la función.

**Extensión con Biopython / PDB (hemoglobina, 1GZX).** La cadena α tiene 141 aa (15 126,2 Da, 46,1 % de residuos hidrofóbicos) y la β, 146 aa (15 867,0 Da, 44,5 %). Los registros `HELIX` y `SHEET` del PDB dan la estructura secundaria: las globinas se pliegan en hélices α, sin apenas láminas β. Con la superficie accesible calculada por `ShrakeRupley`, en la cadena β **57 de 65 residuos hidrofóbicos (88 %) están enterrados**, frente a solo **4 de 29 cargados (14 %)**, lo que confirma que el interior depende de residuos hidrofóbicos. La mutación Glu6Val de la anemia falciforme afecta a un residuo de superficie (según la literatura; el notebook calcula su accesibilidad): no desestabiliza el pliegue, pero crea un parche hidrofóbico expuesto que favorece la agregación de la hemoglobina desoxigenada.

### 3.6. Ejercicio 6. Actividad integradora: del ADN a la proteína

Se tomó TP53 (ENST00000269305.9, 1182 nt) y el pipeline `pipeline_dogma_central` ejecuta tres etapas, informando de cada una con mensajes `[Replicación]`, `[Transcripción]` y `[Traducción]`:

- **Replicación:** muestra las hebras parentales y las dos moléculas hijas, y comprueba que son idénticas a la original.
- **Transcripción:** sintetiza el ARNm leyendo la hebra molde y comprueba que coincide con la codificante con T→U.
- **Traducción:** traduce hasta el codón de paro y avisa si falta, si es prematuro o si el último codón está incompleto.

Antes de empezar valida el CDS (longitud múltiplo de 3, ATG inicial, codón de paro final) y emite un `[AVISO]` si algo no cuadra. El resultado es una proteína de 393 aa, la longitud real de p53 humana.

**Reflexión: ¿qué punto es más vulnerable?** Por frecuencia, la traducción es la que más se equivoca (del orden de un error cada 10 000 codones), seguida de la transcripción, y la replicación la que menos (del orden de un error cada mil millones de bases tras la corrección). Pero los errores de traducción y transcripción son transitorios: afectan a copias individuales y hay muchas otras correctas. Los de replicación quedan **fijados en el ADN**, se heredan y afectan a todos los ARNm y proteínas del gen. Por eso, para la función de la proteína, el punto más crítico es la replicación (una mutación sin sentido o un cambio de marco en TP53 anula la proteína en todas las copias), y es donde la célula concentra más mecanismos de control (corrección de pruebas y reparación de apareamientos erróneos).

## 4. Conclusión

- Biopython reproduce cada paso del dogma central con pocas funciones, pero no orienta las hebras por sí mismo: `transcribe()` da un ARNm correcto solo si recibe la hebra codificante en 5'→3'.
- Los datos reales confirman la teoría: *lacZ* da una proteína de 1024 aa, TP53 de 393 aa y FGFR2 produce 46 transcritos con proteínas de 51 a 822 aa a partir de un solo gen.
- En la hemoglobina el interior está formado casi solo por residuos hidrofóbicos, lo que explica por qué un cambio interno es más grave que uno de superficie.
- Validar las secuencias (marco de lectura, ATG, codón de paro) evita interpretar como isoformas lo que son anotaciones parciales.
