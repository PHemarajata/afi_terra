# เอกสารอธิบายการวิเคราะห์ข้อมูล 16S rRNA สำหรับทีมห้องปฏิบัติการ AFI

> **วัตถุประสงค์ของเอกสารฉบับนี้**
>
> ทีมห้องปฏิบัติการ AFI เป็นเจ้าของโครงการและรับหน้าที่นำเขียนต้นฉบับ (manuscript) เพื่อตีพิมพ์ ทีม bioinformatics ได้ทำการวิเคราะห์ข้อมูล 16S rRNA และเตรียมร่างต้นฉบับภาคบรรยายผลและวิธีการเอาไว้แล้ว เอกสารฉบับนี้สรุปเป็นภาษาไทยว่า**เราทำอะไร ทำไม และผลที่ได้คืออะไร** เพื่อให้ทีมห้องปฏิบัติการอ่านแล้วเข้าใจทั้งกระบวนการก่อนนำไปเขียนต้นฉบับฉบับสมบูรณ์
>
> เอกสารต้นฉบับภาษาอังกฤษ (ฉบับเต็มและฉบับย่อ) พร้อมตารางและรูปประกอบ อยู่ในไฟล์ต่อไปนี้
>
> - `MANUSCRIPT-FINAL-DRAFT.md` — ต้นฉบับเต็ม (ประมาณ 8,500 คำภาษาอังกฤษ)
> - `MANUSCRIPT-CONDENSED-DRAFT.md` — ต้นฉบับฉบับย่อสำหรับวารสารขนาดสั้น เช่น JCM (ประมาณ 3,900 คำ)
> - `APPENDIX-VALIDATION-PANEL.md`, `APPENDIX-STUDY-SAMPLES.md`, `APPENDICES.xlsx` — ตารางผลรายตัวอย่าง
> - `figure_a_sankey.html` / `figure_a_sankey.png` — รูปที่ 1 แสดงการกระจายของจุลชีพก่อน/หลังการกรอง
>
> เอกสารฉบับนี้เขียนวันที่ 2026-05-12

---

## 1. ปัญหาที่เรากำลังพยายามตอบ

โครงการนี้ศึกษาผู้ป่วย Acute Febrile Illness (AFI) หรือผู้ป่วยที่มาด้วยอาการไข้เฉียบพลัน ซึ่งในภาคตะวันออกเฉียงเหนือของประเทศไทย มีสาเหตุได้หลากหลายมาก ตั้งแต่

- **เชื้อแบคทีเรียที่เพาะขึ้นยากหรือเพาะไม่ขึ้นบนวุ้นเลือดธรรมดา** เช่น *Orientia tsutsugamushi* (สครับไทฟัส), *Rickettsia* (โรคไข้พบเห็บ), *Leptospira* (เลปโตสไปโรซิส), *Brucella*
- **เชื้อแบคทีเรียที่เพาะขึ้นได้แต่ต้องอาศัยอาหารเลี้ยงเชื้อพิเศษหรือเวลาเพาะนานกว่าปกติ** เช่น *Burkholderia pseudomallei* (เมลิออยโดสิส)
- **เชื้อแบคทีเรียที่เพาะขึ้นได้ตามปกติ** เช่น *Escherichia coli*, *Klebsiella*
- **ไวรัส โปรโตซัว และเชื้อราอื่น ๆ** เช่น Dengue, Plasmodium

ในกลุ่มตัวอย่างที่เราวิเคราะห์ ผู้ป่วยทุกคนได้รับการเพาะเชื้อจากเลือด (blood culture) ในระหว่างที่นอนรพ. และพบว่า**อย่างน้อย 1 ขวดเพาะเชื้อให้ผลบวกจากเครื่องตรวจอัตโนมัติ** (ขวดมีสัญญาณว่ามีจุลชีพเติบโต) แต่เมื่อนำมาเพาะต่อบนวุ้นเลือดในห้องปฏิบัติการ (subculture) **ไม่พบเชื้อขึ้น** ในระยะเวลาที่กำหนด (5–7 วันโดยปกติ) และโรงพยาบาลที่เก็บตัวอย่างไม่ได้ทำการเพาะเชื้อในสภาวะไม่มีออกซิเจน (anaerobic culture) ตามโพรโทคอลของหน่วยงาน

**คำถามวิจัยคือ: จุลชีพอะไรที่ทำให้ขวดเพาะเชื้อเป็น positive แต่ไม่ขึ้นในการเพาะต่อ?** การสุ่มลำดับ 16S rRNA gene จากตัวอย่างเลือดของผู้ป่วยอาจช่วยตอบคำถามนี้ได้ เพราะการสุ่มลำดับสามารถตรวจพบ DNA ของจุลชีพได้แม้ว่าจุลชีพนั้นจะเพาะไม่ขึ้น

**ข้อควรระวังที่สำคัญตั้งแต่ต้น:** ตัวอย่าง 16S ที่เราวิเคราะห์ **มาจากเลือดของผู้ป่วยในระหว่างเข้ารับการรักษา ไม่ใช่จากขวดเพาะเชื้อโดยตรง** ดังนั้นจุลชีพที่ตรวจพบใน 16S **อาจ**เป็นสาเหตุของสัญญาณบวกในขวดเพาะเชื้อ แต่ก็อาจไม่ใช่ก็ได้ ความเชื่อมโยงเป็น correlative ไม่ใช่ causal

---

## 2. ข้อมูลที่เราใช้

เรามีตัวอย่างสองชุดที่ต้องวิเคราะห์ต่างกัน

### 2.1 ชุด validation panel (43 ตัวอย่าง + NTC 5 ตัวอย่าง = 48 ตัวอย่างรวม)

เป็นชุดตัวอย่างที่ทราบผลล่วงหน้าว่าควรเจอเชื้ออะไร เพื่อใช้ทดสอบว่า pipeline ของเราตรวจจับเชื้อได้ถูกต้องหรือไม่ ประกอบด้วย

- **33 ตัวอย่างผู้ป่วย** ที่ทราบ diagnosis แล้วจากห้องปฏิบัติการอ้างอิง: *E. coli* 5, *Orientia* 6, *Rickettsia* 4, *Leptospira* 4, *B. pseudomallei* 5, *S. pneumoniae* 3, *S. suis* 3, *Coxiella* 2, *Yersinia* 1
- **10 ตัวอย่าง positive control (PC)**
  - **PC_SINGLE 4 ตัว** (เชื้อเดียวต่อตัวอย่าง): *E. coli*, *P. aeruginosa*, *S. pneumoniae*, *S. suis*
  - **PC_MIX8 5 ตัว** เป็น replicate ของ ZymoBIOMICS Microbial Community Standard ซึ่งประกอบด้วยเชื้อ 8 ชนิด: *Bacillus*, *Enterococcus*, *Escherichia*, *Limosilactobacillus*, *Listeria*, *Pseudomonas*, *Salmonella*, *Staphylococcus*
  - **MIXED4 1 ตัว** ผสมเชื้อ 4 ชนิด: *E. coli*, *P. aeruginosa*, *S. pneumoniae*, *S. suis*
- **5 ตัวอย่าง negative template control (NTC)** เป็นน้ำสะอาดที่ผ่านขั้นตอน DNA extraction และ library prep เหมือนตัวอย่างผู้ป่วย ใช้ตรวจการปนเปื้อนจากน้ำยาหรือชุดทดสอบ

### 2.2 ชุดตัวอย่างศึกษา study cohort (86 ตัวอย่าง)

เป็นตัวอย่างจริงจากผู้ป่วย AFI ที่เพาะเลือดได้ positive bottle แต่ไม่ขึ้นเชื้อใน subculture (กลุ่มเป้าหมายของงานวิจัย) จัดอยู่ใน 5 run การ sequencing (1_and_2, 3, 4_and_5, 6_and_7, 8_and_9)

> **หมายเหตุสำคัญ:** ชุดศึกษามี **86 ตัวอย่าง** ที่ผ่าน pipeline (มีไฟล์ผลลัพธ์) ในจำนวนนี้ 15 ตัวอย่างไม่พบจุลชีพใด ๆ เลย (zero detection — น่าจะเกิดจาก DNA extraction / library prep ล้มเหลว) ส่วน 71 ตัวอย่างพบจุลชีพอย่างน้อย 1 ชนิด

---

## 3. ภาพรวม pipeline ขั้นตอนการวิเคราะห์

ขั้นตอนการประมวลผลตั้งแต่ raw reads จนถึงรายงานผลทำผ่าน WDL workflow บน Terra/Cromwell (ซอร์สโค้ดอยู่ที่ https://github.com/PHemarajata/afi_terra) สรุปเป็นภาพรวมได้ดังนี้

```
ตัวอย่างเลือด → DNA extraction → ตัด primer/adapter → ลำดับ 16S V1-V3
                                                            ↓
  raw reads (paired-end FASTQ)
                                                            ↓
        [ขั้นที่ 1] NCBI Scrubber — ลบลำดับมนุษย์ออก
        [ขั้นที่ 1.5] fastp — กรองคุณภาพและตัด adapter
                                                            ↓
        [ขั้นที่ 2] Centrifuger v1.1.0 — จัดประเภทเชื้อระดับ genus
              → ผลลัพธ์: kreport ที่บอกว่ามี genus อะไรกี่ reads
                                                            ↓
        [ขั้นที่ 3] minimap2 — alignment เฉพาะกลุ่ม Rickettsiales
              → ใช้ reference panel 7 จีโนม (Orientia 2 strain, Rickettsia 3
                สายพันธุ์, Anaplasma 1, Ehrlichia 1)
              → ผลลัพธ์: mapped_reads และ breadth of coverage
                                                            ↓
        [ขั้นที่ 4] รวมผลและเปรียบเทียบกับ NTC ของ run นั้น (NCmax)
              → กำหนด call: Detected / Confirmed / Probable / Negative
                                                            ↓
        [ขั้นที่ 5] V4 Decontamination Filter — กรองจุลชีพปนเปื้อน
                                                            ↓
        ผลรายงานสุดท้าย: เก็บไว้ในไฟล์ .calls.tsv และ .taxa_evidence.tsv
```

### 3.1 ขั้นที่ 1: ลบลำดับมนุษย์ (NCBI Scrubber + fastp)

เนื่องจากเลือดของผู้ป่วยมี DNA ของมนุษย์ปริมาณมาก จึงต้องลบออกก่อน ไม่งั้นการวิเคราะห์ขั้นต่อไปจะช้าและไม่มีประสิทธิภาพ ใช้ NCBI SRA Human Scrubber (Katz et al., 2021) ตามด้วย fastp (Chen et al., 2018) เพื่อตัด adapter และกรอง read คุณภาพต่ำ

### 3.2 ขั้นที่ 2: Centrifuger — จัดประเภทเชื้อ

Centrifuger เป็นโปรแกรมที่นำลำดับ DNA ที่เหลือไปเทียบกับฐานข้อมูลจีโนมของแบคทีเรียทั้งหมด แล้วบอกว่าแต่ละ read น่าจะเป็นของเชื้ออะไร ผลลัพธ์รวมออกมาเป็น **kreport** ซึ่งเป็นตารางสรุปว่าตัวอย่างนี้มี genus อะไรกี่ reads และ species อะไรกี่ reads

**เกณฑ์การเรียก "Detected" สำหรับ genus:**
- ต้องมี read ≥ 500 reads
- AND ต้องมีจำนวน read ≥ 5 เท่าของค่า NCmax ของ genus นั้นใน run เดียวกัน (เพื่อตัดสิ่งปนเปื้อนจากน้ำยา)

ทั้งสองเงื่อนไขต้องเป็นจริงพร้อมกัน

> **คำว่า "Centrifuger" vs "Centrifuge"**: ทั้งสองเป็นซอฟต์แวร์คนละตัว Centrifuger เป็นรุ่นใหม่กว่า (Song & Langmead 2024) ที่พัฒนาต่อจาก Centrifuge (Kim et al. 2016) — pipeline ของเราใช้ Centrifuger ดังนั้นไฟล์ผลลัพธ์มีชื่อว่า `*.centrifuger.kreport.tsv` (ลงท้ายด้วย -er) ในเอกสารต้นฉบับภาษาอังกฤษเราอ้างอิงทั้งสองตัว

### 3.3 ขั้นที่ 3: minimap2 alignment สำหรับ Rickettsiales

**ทำไมต้องมีขั้นนี้?** เพราะการตรวจจับเชื้อกลุ่ม Rickettsiales (Orientia + Rickettsia) ที่ระดับ V1-V3 ของ 16S rRNA ทำได้ยาก ลำดับ DNA ระหว่าง Orientia กับ Rickettsia ในบริเวณนี้คล้ายกันมาก ทำให้ Centrifuger บางครั้งแยกไม่ออกหรือเรียกเป็น "low-confidence call" Rickettsiales เป็น**สาเหตุของ AFI ที่พลาดไม่ได้** (cannot-miss diagnosis) ในประเทศไทยภาคตะวันออกเฉียงเหนือ จึงต้องเพิ่มขั้นตอน rescue โดยเฉพาะ

minimap2 จะนำ read ทั้งหมดไปจัด alignment กับ reference panel ของ Rickettsiales (7 จีโนม) แล้วคำนวณ
- **mapped reads** = จำนวน read ที่ align ได้
- **breadth of coverage** = สัดส่วนของ reference genome ที่ถูก cover ด้วย read (เช่น 0.30 = cover 30% ของ reference)

จากนั้นใช้เกณฑ์สองชั้น (two-tier framework) ในการเรียกผล

#### Tier 1 — Confirmed (ระบุ genus ได้ชัดเจน)
- mapped reads ≥ 100
- AND breadth ≥ 0.25
- AND mapped reads ≥ 5 เท่าของ NCmax สำหรับ alignment
- รายงานเป็น *Orientia* หรือ *Rickettsia* โดยเฉพาะ

#### Tier 2 — Probable ("Rickettsiales detected" ระดับ order)
- mapped reads ≥ 50
- AND breadth ≥ 0.20
- AND mapped reads > NCmax (ไม่ต้องถึง 5 เท่า)
- รายงานเป็น **"Rickettsiales detected (genus uncertain; recommend confirmatory species-specific qPCR)"**
- มีนัยทางคลินิกคือ ผู้ป่วยอาจมีการติดเชื้อในกลุ่ม Rickettsiales ซึ่งตอบสนองต่อ doxycycline เป็นยาตัวเลือกแรกอยู่แล้วในผู้ป่วย AFI ในประเทศไทย

#### ไม่ผ่านเกณฑ์ — Negative
ผลถูกเก็บไว้ในไฟล์เพื่อความโปร่งใส (สามารถตรวจสอบย้อนหลังได้ว่าทำไมไม่ผ่าน) แต่ไม่รายงานทางคลินิก

### 3.4 ขั้นที่ 4: NTC subtraction (per-run NCmax)

ในการวิเคราะห์ตัวอย่างชีวภาพมวลต่ำ (low-biomass — เลือดถือเป็นตัวอย่างชนิดนี้) ปัญหาที่หลีกเลี่ยงไม่ได้คือ**การปนเปื้อนจากน้ำยา ชุดสกัด DNA และอุปกรณ์ห้องปฏิบัติการ** จุลชีพหลายชนิดที่ดูเหมือนจะตรวจพบจริง ๆ แล้วเป็นการปนเปื้อน

วิธีจัดการคือใส่ NTC (น้ำเปล่าผ่านขั้นตอน extraction และ library prep เหมือนตัวอย่างผู้ป่วย) เข้าไปใน run และใช้ผลของ NTC เป็น "พื้นหลัง" (background) สำหรับเปรียบเทียบ

**กฎการคำนวณ NCmax (จาก `scripts/build_ntc_background.py`):**
- สำหรับแต่ละ genus และแต่ละ run คำนวณ **ค่า read สูงสุด** (maximum) ที่พบใน NTC ของ run นั้น
- คำนวณแยกระหว่าง Centrifuger reads (`cfr_ntc_reads`) และ Minimap2 alignment reads (`align_ntc_reads`)
- ตัวอย่างผู้ป่วยใน run นั้นจะต้องมี read มากกว่า NCmax คูณ fold ที่กำหนด (5×) จึงจะนับว่า "Detected"

> **ข้อสังเกตที่พบระหว่างการทบทวน:** ใน run 6_and_7 พบว่า NTC2_ExDw_S13_L001 มีการปนเปื้อนสูงมาก (Leptospira 78,691 reads, Burkholderia 106,204 reads, Brevundimonas 188,778 reads) อย่างไรก็ตามตัวอย่างผู้ป่วยใน run เดียวกันใน `.calls.tsv` แสดง `ntc_reads = 0` สำหรับ Leptospira ซึ่งหมายความว่าระบบ assignment ของ run-ID อาจจะแยก NTC ตัวที่ปนเปื้อนนี้ไปอยู่ใน NTC pool คนละกลุ่มกับตัวอย่างที่ได้รับผลกระทบ การจัดกลุ่ม NTC ใน pipeline ขึ้นอยู่กับ Terra-level run-ID ไม่ใช่โฟลเดอร์ที่จัดไว้ ในต้นฉบับ section §2.4.4 มีอธิบายรายละเอียดและแนะนำให้ระบุชัดเจนใน deployment SOP

### 3.5 ขั้นที่ 5: V4 Decontamination Filter

หลังจากผ่าน pipeline ขั้นที่ 1–4 แล้ว เรายังมี "Detected" calls อีกจำนวนหนึ่งที่น่าจะเป็นการปนเปื้อน เพราะ NTC subtraction ไม่สามารถจัดการกับการปนเปื้อนทั้งหมดได้ (โดยเฉพาะการปนเปื้อนที่ NTC ตัวที่ใช้คำนวณ NCmax ก็มีปนเปื้อนเหมือนกัน ทำให้ค่า NCmax ต่ำกว่าที่ควรจะเป็น) เราจึงเพิ่ม **decontamination filter** อีกชั้นหนึ่ง

V4 filter ประกอบด้วย 4 tier และมีข้อยกเว้นพิเศษอีก 1 ข้อ

#### Tier A — จุลชีพปนเปื้อนความเชื่อมั่นสูง (high-confidence kit/skin/water contaminants)

ตัด genus ต่อไปนี้ออกทุกครั้งที่พบในตัวอย่างผู้ป่วยหรือศึกษา (รวม 11 genus):

*Pseudomonas, Ralstonia, Bradyrhizobium, Sphingomonas, Stenotrophomonas, Methylobacterium, Acinetobacter, Cutibacterium, Staphylococcus, Corynebacterium, Brevundimonas*

genus เหล่านี้ปรากฏใน landmark studies เกี่ยวกับการปนเปื้อนในตัวอย่าง low-biomass อย่างน้อย 5 papers (Salter 2014, Glassing 2016, Lauder 2016, de Goffau 2018, Tan 2023) และในข้อมูลของเราเองก็พบว่ามักจะตรวจพบใน NTC ในปริมาณมาก

> **ข้อยกเว้น "PC bypass":** สำหรับตัวอย่าง positive control ที่จุลชีพในกลุ่ม Tier A เป็น "expected spike-in" จะไม่ถูกตัด เช่น *Pseudomonas* ตัดออกในตัวอย่างผู้ป่วย แต่ใน P-aeru_S5_L001 และ PC_MIX8 *Pseudomonas* คือเชื้อที่เราจงใจใส่ลงไป จึงต้องคงเอาไว้ (อ่านรายละเอียดใน section §4 ของเอกสารฉบับนี้)

#### Tier A exception — Burkholderia species-level safeguard (สำคัญมาก)

***Burkholderia* เป็น kit contaminant ในกลุ่ม Tier A ทั่วไป แต่ภายในสกุลนี้มี *B. pseudomallei* ซึ่งเป็นสาเหตุของเมลิออยโดสิส ซึ่ง endemic ในประเทศไทยและเป็น cannot-miss diagnosis** เราจึงไม่สามารถตัด *Burkholderia* ทิ้งทั้งหมด

วิธีจัดการคือ **เปิดไฟล์ kreport ของ Centrifuger อ่านระดับ species** (rank `S`) แล้วดูว่ามีกี่ reads ที่ถูก assign ให้เป็น *Burkholderia pseudomallei* โดยเฉพาะ (NCBI taxonomy 28450)

- ถ้า species-level reads ≥ 500 AND species-level reads > NCmax species-level ของ run → **เก็บไว้** (preserved) และรายงานเป็น *B. pseudomallei* species-confirmed
- ถ้าไม่ครบเงื่อนไข → **ตัดทิ้ง** ทั้งหมด (ถือว่า genus signal เป็น contamination จาก *B. cepacia* complex หรือ species อื่น)

ใน validation panel: ตัวอย่างที่คาดว่าเป็น *B. pseudomallei* 3 ตัว (09502813_S2_L001, 09-0-02165, 09700912_S3_L001) มี species-level reads 13,744 / 20,061 / 65,016 ซึ่งผ่านเกณฑ์ทั้งหมด

ใน study cohort: ไม่มีตัวอย่างใดผ่านเกณฑ์ (รายละเอียดใน section 5.2 ของเอกสารนี้)

#### Tier B — NTC-only organisms

9 genus ที่พบเฉพาะใน NTC และไม่เคยพบในตัวอย่างผู้ป่วยในข้อมูลของเราเลย: *Cereibacter, Thioclava, Bdellovibrio, Saltatorellus, Pseudogemmobacter, Minisyncoccus, Rhodoluna, Microbacterium, Arcanobacterium*

#### Tier 1 — ultra-low abundance noise

16 genus ที่มี median abundance < 0.5% และไม่น่าจะเป็นเชื้อก่อโรค: *Shigella, Metapseudomonas, Stutzerimonas, Capsulimonas, Chamaesiphon, Chloroflexus, Flavihumibacter, Hymenobacter, Limnoglobus, Methylovirgula, Microvirga, Pelagovum, Pseudonocardia, Rufibacter, Salmonella, Spirosoma*

(*Mycoplasmopsis* และ *Nitrospira* ไม่อยู่ในรายการนี้ภายใต้ V4 เนื่องจากข้อมูลจริงแสดงปริมาณสูงกว่าเกณฑ์ ultra-low-abundance)

#### Tier 2 — marginal organisms

27 genus ที่มี median abundance 0.5–2.0% เก็บไว้เฉพาะกรณีที่พบใน ≥ 2 ตัวอย่าง AND แต่ละ detection มี abundance ≥ 1.0%

#### PC spike-in bypass (สำคัญสำหรับ validation panel)

สำหรับตัวอย่าง positive control ระบบจะตรวจสอบรายการ "expected spike-in organisms" ของตัวอย่างนั้น ๆ และ**ข้าม Tier A removal** สำหรับ genus ที่อยู่ในรายการ ตารางต่อไปนี้แสดง expected organisms ของแต่ละ PC

| PC sample | ประเภท | จุลชีพที่คาดว่ามี (genera) |
|---|---|---|
| E-coli_S4_L001 | PC_SINGLE | *Escherichia* |
| P-aeru_S5_L001 | PC_SINGLE | *Pseudomonas* |
| S-pneumo_S2_L001 | PC_SINGLE | *Streptococcus* |
| S-suis_S3_L001 | PC_SINGLE | *Streptococcus* |
| PC_MIX8 (5 replicates) | PC_MIX8 | *Bacillus, Enterococcus, Escherichia, Limosilactobacillus, Listeria, Pseudomonas, Salmonella, Staphylococcus* (ZymoBIOMICS Standard) |
| Mixed_S6_L001 | MIXED4 | *Escherichia, Pseudomonas, Streptococcus* |

ถ้าไม่มี PC bypass ตัวอย่าง P-aeru_S5_L001 จะ "fail" เพราะ *Pseudomonas* (เชื้อที่ใส่ลงไปจงใจ) ถูกตัดทิ้งโดย Tier A เมื่อมี bypass แล้ว PC ทั้ง 10 ตัวอย่างผ่านเกณฑ์ 100% และผลรวมของ sample-level analytical performance ขึ้นไปอยู่ที่ 31/43 = 72.1% ซึ่งตรงกับตัวเลขที่รายงานใน APHL validation report เดิม

---

## 4. ผลการ validate ของ pipeline (validation panel)

นี่คือผลการทดสอบ pipeline กับ validation panel 48 ตัวอย่าง

### 4.1 ผลรวม (overall accuracy)

| หมวด | Concordant | จำนวนทั้งหมด | อัตรา | 95% CI |
|---|---|---|---|---|
| Clinical (เจอ target organism และผ่าน V4) | 21 | 33 | **63.6%** | 46.6–77.8% |
| Positive controls (เจอ spike-in ครบทุก organism) | 10 | 10 | **100%** | 72.2–100% |
| NTCs (ไม่มี TAC bacterial target ใด ๆ ผ่าน V4) | 5 | 5 | **100%** | 56.6–100% |
| **Sample-level analytical performance (clinical + PC)** | **31** | **43** | **72.1%** | 57.3–83.3% |
| **Overall validation accuracy (ทุกหมวด)** | **36** | **48** | **75.0%** | 61.2–85.1% |

ตัวเลข **72.1% sample-level analytical performance** คือเลขที่จะใช้รายงานในต้นฉบับ ตรงกับ APHL bioinformatic validation report เดิม

### 4.2 ผลแยกตาม organism

| Expected organism | n | Concordant | Sensitivity | สาเหตุของ discordant |
|---|---|---|---|---|
| *E. coli* | 5 | 5 | 100% | — |
| *O. tsutsugamushi* | 6 | 6 | 100% | (rescued ผ่าน Tier 1/2) |
| *Rickettsia* spp | 4 | 2 | 50% | 2 ตัว pre-sequencing failure (00126_S6, 00369_S1: zero taxa) |
| *Leptospira* spp | 4 | 2 | 50% | 2 ตัว abundance-driven (organism อื่นครอบครอง signal) |
| *B. pseudomallei* | 5 | 3 | 60% | 2 ตัว abundance-driven |
| *S. pneumoniae* | 3 | 1 | 33% | 2 ตัว abundance-driven; V1-V3 แยก S. pneumoniae vs S. suis ไม่ได้ |
| *S. suis* | 3 | 1 | 33% | เหมือนข้างต้น |
| *C. burnetii* | 2 | 1 | 50% | 1 ตัว pre-sequencing failure (25800370_S9: zero taxa) |
| *Yersinia* spp | 1 | 0 | 0% | n=1 ไม่เพียงพอสำหรับสรุป |

### 4.3 ข้อสังเกตที่สำคัญจาก validation panel

1. **3 ตัวอย่าง pre-sequencing failure** (00126_S6, 00369_S1, 25800370_S9) มี zero taxa — ปัญหาน่าจะอยู่ที่ DNA extraction / library prep / sequencing depth ไม่ใช่ที่ classifier
2. **6 ตัวอย่างของ Streptococcus** (S. pneumoniae 3 + S. suis 3) — V1-V3 ไม่สามารถแยก species ได้ คงต้องรายงานในระดับ genus เท่านั้น
3. **00618_S7_L001 (Orientia expected)** — `.calls.tsv` แสดง `call = Probable` (Tier 2 order-level rescue triggered) นับเป็น concordant ตามกรอบ two-tier framework
4. **V4 filter ไม่ตัด target organism ใด ๆ ใน validation panel** — เพราะตัวอย่าง validation มักจะมี expected organism เป็น dominant signal (abundance สูง) อยู่แล้ว ไม่มี Tier A contaminant ที่ตรวจพบในระดับที่มีผลต่อ target

### 4.4 ผลของ control validity

- **PC_MIX8 5 ตัว**: ทุกตัวตรวจพบทั้ง 8 expected organisms ครบ (40/40 detection) เมื่อใช้ PC bypass
- **PC_SINGLE 4 ตัว**: ทุกตัวตรวจพบ expected organism ครบ
- **MIXED4 1 ตัว**: ตรวจพบทั้ง 3 expected genera ครบ
- **NTC 5 ตัว**: ไม่มี TAC bacterial target ใด ๆ ผ่าน V4 (specificity 100%)
- **Inter-run reproducibility 100%** ทั้ง 9 sequencing run

---

## 5. ผลการวิเคราะห์ study cohort

นี่คือผลการ apply pipeline ที่ผ่านการ validate แล้วกับตัวอย่างศึกษาจริง 86 ตัวอย่าง

### 5.1 ผลรวม pre-sequencing failures และจำนวนตัวอย่างที่มีผล

- **71 / 86 ตัวอย่าง (82.6%)** มี ≥ 1 positive call (Detected, Confirmed, หรือ Probable)
- **15 / 86 ตัวอย่าง (17.4%)** มี zero detection — เป็น pre-sequencing failure แนะนำให้พิจารณาทำซ้ำหรือใช้ diagnostic modality อื่น

### 5.2 ข้อค้นพบสำคัญ

มี 3 ข้อค้นพบที่ทีมต้องเข้าใจให้ชัดก่อนเขียนต้นฉบับ

#### ข้อค้นพบที่ 1 (สำคัญที่สุด): พบ Rickettsiales rescue ใน 11 / 86 ตัวอย่าง (12.8%)

หลักฐานนี้มาจากแถว Minimap2 alignment rescue ใน `.calls.tsv` (rows ที่มี `source = alignment` และ `call ∈ {Confirmed, Probable, Detected}`)

11 ตัวอย่างที่มี Rickettsiales rescue:

| ตัวอย่าง | Run | Genus | mapped reads | breadth | Tier | NTC reads | Confidence |
|---|---|---|---|---|---|---|---|
| 16901195_S5_L001 | 8_and_9 | Orientia | 14,816 | 0.3235 | **Tier 1 (Confirmed)** | 1,288 | HIGH |
| 22600400_S6_L001 | 8_and_9 | Orientia | 9,747 | 0.2181 | Tier 2 | 1,288 | HIGH |
| 23900356_S5_L001 | 4_and_5 | Rickettsia | 441,173 | 0.3046 | Tier 2 | 189,002 | MODERATE |
| 23200430_S6_L001 | 4_and_5 | Rickettsia | 94,142 | 0.2189 | Tier 2 | 53,908 | LOW |
| 23200519_S12_L001 | 6_and_7 | Orientia | 2,492 | 0.2175 | Tier 2 | 1,143 | MODERATE |
| 09801652_S5_L001 | 6_and_7 | Orientia | 3,727 | 0.2282 | Tier 2 | 1,143 | MODERATE |
| 25800718_S8_L001 | 6_and_7 | Orientia | 4,143 | 0.2248 | Tier 2 | 2,202 | LOW |
| 10100409_S9_L001 | 8_and_9 | Orientia | 2,385 | 0.2221 | Tier 2 | 1,288 | LOW |
| 16601093_S7_L001 | 8_and_9 | Orientia | 2,284 | 0.2201 | Tier 2 | 1,288 | LOW |
| 23200736_S8_L001 | 8_and_9 | Orientia | 1,981 | 0.2141 | Tier 2 | 1,288 | LOW |
| 23900752_S9_L001 | 8_and_9 | Orientia | 1,684 | 0.2428 | Tier 2 | 1,288 | LOW |

**ความหมายทางคลินิก:** Rickettsiales เป็นเชื้อ obligate intracellular ที่**เพาะไม่ขึ้นบนวุ้นเลือดธรรมดา** การที่ตรวจพบ DNA ของ Rickettsiales ใน 11/86 ตัวอย่าง (ประมาณ 13% ของ cohort) สอดคล้องกับ epidemiology ของ AFI ใน NE Thailand เป็นอย่างยิ่ง — Rickettsiales (โดยเฉพาะ scrub typhus) เป็นสาเหตุของ AFI ที่พบบ่อยที่สุดในพื้นที่นี้ และเป็นคำตอบที่ขับเคลื่อนคลินิก (doxycycline empirical therapy)

**Confidence ปรับตามอัตราส่วนระหว่าง sample reads กับ NTC reads:**
- HIGH = sample มี read สูงกว่า NTC มาก (มาก ≥10 เท่า) หรือเป็น Tier 1 confirmed
- MODERATE = sample > NTC แต่ headroom ไม่มาก
- LOW = sample ใกล้ NTC (NTC carries comparable signal) — ควรยืนยันด้วย qPCR

**คำแนะนำ:** ทั้ง 11 ตัวอย่างควรทำการยืนยันด้วย species-specific qPCR (Orientia-tsutsugamushi-specific หรือ Rickettsia 17 kDa antigen gene targets) ก่อนรายงานผลทางคลินิก

#### ข้อค้นพบที่ 2: ไม่พบ B. pseudomallei ในระดับ species ใด ๆ ใน study cohort

V4 filter parse Centrifuger kreport ที่ระดับ species (rank `S`) สำหรับ *Burkholderia pseudomallei* และจะเก็บไว้เป็น *B. pseudomallei* เฉพาะเมื่อ species-level reads ≥ 500 และสูงกว่า NTC ในรอบเดียวกัน

ในตัวอย่าง 23200430_S6_L001 ซึ่งมี *Burkholderia* genus-level reads = 54,126 เมื่อ parse ที่ระดับ species พบว่า

| Burkholderia species | reads |
|---|---|
| *B. contaminans* | 3,266 |
| *B. cenocepacia* | 2,988 |
| *B. multivorans* | 2,170 |
| *B. sola* | 1,443 |
| *B. cepacia* | 1,335 |
| (B. cepacia complex อื่น ๆ ) | ~ 2,000 รวม |
| ***Burkholderia pseudomallei*** | **8** |
| *B. thailandensis* | 5 |
| *B. mallei* | 3 |

นั่นคือ genus-level signal ของ *Burkholderia* ในตัวอย่างนี้ ส่วนใหญ่เป็น ***B. cepacia* complex** ซึ่งเป็น kit/water contaminant ที่ทราบกันดี ส่วน *B. pseudomallei* species เพียง **8 reads** เท่านั้น ซึ่งน้อยกว่า detection threshold (500) และน้อยกว่า NTC ใน run เดียวกัน (10 reads)

V4 filter ที่แก้ไขใหม่จึงตัด *Burkholderia* ในตัวอย่างนี้ทิ้ง (เพราะไม่ผ่าน species-level safeguard) และ**ไม่มีตัวอย่าง study ใดผ่านเกณฑ์ *B. pseudomallei* species-level** ตัวอย่าง *B. pseudomallei* ใน validation panel (3 ตัว 13,744–65,016 species reads) ยังคงเป็น "ของจริง"

**ความหมายทางคลินิก:** Cohort นี้ไม่มี melioidosis ที่ detect ได้ด้วย 16S V1-V3 — แต่ assay นี้ตรวจ *B. pseudomallei* ได้ดี (แสดงใน validation) เพียงแต่ในกลุ่มผู้ป่วยที่เราศึกษา ไม่มีผู้ป่วย melioidosis ในจำนวนนี้

#### ข้อค้นพบที่ 3: Mycoplasmopsis เป็น candidate signal ของเชื้อ fastidious ที่ใหญ่ที่สุดใน cohort

***Mycoplasmopsis*** ตรวจพบใน 4 ตัวอย่าง โดยมี read count 537–6,568 reads (mean abundance 39.78%, max 70.55%) V4 filter คงสัญญาณนี้ไว้ (อยู่นอก Tier 1 ultra-low-abundance)

**ความหมายทางคลินิก:** Mycoplasma-class organisms เป็นเชื้อที่
- **ไม่มี cell wall** (no peptidoglycan)
- **fastidious** ต้องการ media ที่ supplemented ด้วย sterol เช่น Mycoplasma broth
- **ช้า** ใช้เวลาประมาณ 1–3 สัปดาห์เพื่อให้เห็น colony

นี่คือ profile ของเชื้อที่อธิบาย "positive bottle / no subculture growth" ได้ดีที่สุดในเชิงชีววิทยา — ขวดเพาะเชื้อมี signal เพราะเชื้อมี metabolic activity ในของเหลว แต่เพาะต่อบนวุ้นเลือดธรรมดาในเวลา 5–7 วัน เชื้อนี้ไม่ขึ้น

**คำแนะนำ:** ทั้ง 4 ตัวอย่างควรยืนยันด้วย Mycoplasma-specific PCR และพิจารณาเพาะเชื้อบน specialized media (sterol-enriched)

### 5.3 Candidate signals อื่น ๆ ที่ตรวจพบ (ต้องยืนยัน)

- **Leptospira (1 ตัวอย่าง):** 09801652_S5_L001 มี 7,294 reads (22.78% ของ sample) แต่**ใน run เดียวกัน NTC2_ExDw_S13_L001 มี 78,691 reads ของ Leptospira (10 เท่าของ sample)** เป็น caveat สำคัญ ควรยืนยันด้วย Leptospira-specific qPCR และ paired serology ก่อนแปลผลใด ๆ
- **Brucella (3 ตัวอย่าง):** mean 0.98%, max 1.96% — อยู่ที่ noise floor ควรถือเป็น candidate รอการยืนยัน serology + qPCR
- **Streptococcus (5 ตัวอย่าง):** mean 12.75% — V1-V3 ไม่แยก S. pneumoniae / S. suis และไม่แยก fastidious vs non-fastidious species ได้

### 5.4 จุลชีพอื่น ๆ ที่ตรวจพบใน cohort (full list อยู่ใน APPENDIX-STUDY-SAMPLES.md)

- *Thermomicrobium* (11 ตัวอย่าง, mean 9.13%) — environmental thermophile, น่าจะเป็น noise
- *Escherichia* (6), *Klebsiella* (4), *Enterobacter* (3) — เชื้อก่อโรค AFI ที่เพาะขึ้นได้ตามปกติ
- *Porphyromonas* (1), *Desulfovibrio* (1) — anaerobic genera ที่ตรวจพบ (แสดงว่า "no anaerobes" ของเอกสารเก่าไม่ถูกต้อง)

---

## 6. การวิเคราะห์เพิ่มเติมที่จัดทำไว้สำหรับทีม

### 6.1 ตารางผลรายตัวอย่าง (appendix files)

ทุก detection ในทุกตัวอย่างมีการบันทึกแบบละเอียดในไฟล์ต่อไปนี้

- **`APPENDIX-VALIDATION-PANEL.md`** (markdown): ตาราง 236 แถว แสดง 20 คอลัมน์ต่อแถว — sample, run, category, expected, total reads, biomass kept, biomass removed, % biomass kept, detected genus, source, reads, % of sample, pipeline NTC reads, Minimap2 rescue info, rescue tier, V4 filter outcome, per-row concordance, failure mode, Final concordance (TAC + V4), notes พร้อมตาราง concordance summary รายตัวอย่างท้ายไฟล์
- **`APPENDIX-STUDY-SAMPLES.md`** (markdown): ตาราง 336 แถว 22 คอลัมน์ เพิ่ม same-run NTC max cross-check, Burkholderia species evidence, confidence, hypothesis class, recommended follow-up
- **`APPENDICES.xlsx`**: เนื้อหาเดียวกันในรูป Excel 4 sheets เพื่อให้ทีมเปิดดูใน spreadsheet ได้สะดวก
- **`DECONTAMINATION-FILTER-REPORT-V4.txt`**: output ดิบจาก V4 filter พร้อม summary statistics และ Burkholderia species evidence รายตัวอย่าง

### 6.2 รูปประกอบ

- **`figure_a_sankey.html`** (interactive) และ **`figure_a_sankey.png`** (รูปนิ่ง) — Sankey diagram แสดงการกระจายของจุลชีพ**ก่อน vs หลัง V4 filter** ใน study cohort สามารถ hover ใน HTML เพื่อดู read count ของแต่ละ flow

### 6.3 ซอร์สโค้ดที่ใช้

- `afi_decontamination_filter_v4.py`: V4 filter ที่ใช้รัน
- `generate_appendices.py`: script สำหรับสร้างตาราง appendix (สามารถรันใหม่ได้หาก input data เปลี่ยน)
- `figure_a_sankey.py`: script สำหรับสร้าง Sankey figure

---

## 7. คำแนะนำเรื่องการทดสอบยืนยัน (confirmatory testing)

จุลชีพที่ตรวจพบในการศึกษานี้เป็น **candidate detections** ทั้งหมด ก่อนรายงานผลทางคลินิกควรยืนยันด้วยวิธีอื่นในแต่ละกลุ่ม

| ตัวอย่าง | จำนวน | คำแนะนำการยืนยัน |
|---|---|---|
| Rickettsiales rescue (Orientia/Rickettsia) | 11 | Species-specific qPCR + paired serology; doxycycline empirical coverage ในผู้ป่วย AFI ในพื้นที่นี้เป็น standard อยู่แล้ว |
| *Mycoplasmopsis* candidate | 4 | Mycoplasma-specific PCR; เพาะเชื้อบน Mycoplasma broth / SP4 medium ถ้ามีตัวอย่างเหลือ |
| *Leptospira* candidate (มี NTC caveat) | 1 | Leptospira-specific qPCR + paired serology; ตรวจสอบ NTC2_ExDw_S13_L001 ใน run 6_and_7 |
| *Brucella* candidate (near noise floor) | 3 | Serology IgM/IgG + Brucella-specific qPCR |
| Pre-sequencing failure | 15 | Repeat specimen collection ถ้าทำได้; พิจารณา alternative diagnostic modality |

---

## 8. ข้อจำกัดของการศึกษาที่ทีมต้องตระหนัก

1. **ตัวอย่าง 16S = เลือดผู้ป่วย ไม่ใช่ขวดเพาะเชื้อ** — เราตรวจจาก patient blood ไม่ใช่จากของเหลวใน blood culture bottle ดังนั้นการเชื่อมโยงระหว่าง 16S detection กับ positive bottle signal เป็น correlative ไม่ใช่ causal
2. **ไม่มีการประเมิน viability** — 16S ตรวจ DNA ไม่ว่าเชื้อจะมีชีวิตอยู่หรือไม่ (รวมถึง dead cells และ VBNC state)
3. **V1-V3 primer bias** — 27F primer set ที่ใช้มีปัญหาในการตรวจจับ anaerobes บางกลุ่ม และ Bifidobacterium / Gardnerella บางสายพันธุ์ ผลลบสำหรับ organism กลุ่มเหล่านี้ไม่ได้แปลว่าไม่มีจริง
4. **Genus-level resolution** — เฉพาะ *Burkholderia* ที่ได้ species-level safeguard organism อื่น ๆ ระบุเฉพาะระดับ genus (เช่น Streptococcus, Leptospira, Brucella)
5. **NTC contamination ใน run บางตัว** — NTC2_ExDw_S13_L001 ใน run 6_and_7 มีปนเปื้อนมาก pipeline จัดการผ่าน run-ID-based NTC pool assignment แต่ควรระบุชัดเจนใน deployment SOP
6. **Sample size 86 ตัวอย่างเป็น pilot** — single-sample (เช่น Leptospira) หรือ few-sample (เช่น Brucella) findings ไม่สามารถสรุป cohort-level prevalence ได้ ควรถือเป็น hypothesis-generating สำหรับ follow-up study
7. **ไม่มี clinical outcome correlation** ในเอกสารปัจจุบัน — ข้อมูล demographics, treatment, response, outcome ของผู้ป่วยยังไม่ได้ link กับผล 16S (ส่วนนี้ทีมต้องเพิ่ม)

---

## 9. สิ่งที่ต้องเพิ่มจากทีมก่อนส่งต้นฉบับตีพิมพ์

ตามตารางท้าย `MANUSCRIPT-FINAL-DRAFT.md` มี placeholders ที่ต้องการ team input อยู่ 18 จุด สรุปประเด็นหลัก:

| ส่วน | ข้อมูลที่ต้องเพิ่ม |
|---|---|
| Title page | รายชื่อผู้แต่ง, affiliation, corresponding author, running title, CoI, funding, author contributions |
| §2.1 Study design | IRB approval, ethics body, ระยะเวลา enrollment, รพ./อำเภอ, inclusion / exclusion criteria |
| §2.2 Sample populations | วิธีการ diagnostic อ้างอิงของแต่ละ pathogen ใน validation panel (PCR? culture? serology?) |
| §2.3 Wet-lab | ปริมาณเลือดที่เก็บ, container, DNA extraction kit, 16S primer sequences, library prep kit, sequencer model, run configuration |
| §2.4.2 | ระบุ Centrifuger reference database identifier (ชื่อ, build date, source URL หรือ accession set) |
| §2.7, §3.2.2, §3.3 | reference ของ deployment SOP สำหรับ per-run QC pass criteria |
| §2.8 Data availability | location ที่ archived ของ raw `.calls.tsv` และ kreport (Google Drive AFI_compare? หรือ data repository accession?) |
| §4.5 Clinical follow-up | ผลการ confirmatory testing (ถ้ามี) สำหรับ candidate detections; clinical outcomes ของผู้ป่วย |
| §7 References | AFI epidemiology references สำหรับ NE Thailand; clinical references สำหรับ Rickettsiales / melioidosis / leptospirosis / brucellosis; CDC TaqMan Array Card AFI panel validation reports |

---

## 10. สรุปสำหรับทีม

โดยสรุป pipeline 16S V1-V3 ของเรามีคุณสมบัติดังนี้

- **เกณฑ์การตรวจจับชัดเจน** — Centrifuger ≥500 reads + 5× NCmax; Minimap2 Tier 1 confirmed ≥100 reads + breadth ≥0.25 + 5× NCmax; Tier 2 probable ≥50 reads + breadth ≥0.20 + > NCmax
- **มี two-tier Rickettsiales framework** — รองรับว่า V1-V3 แยก Orientia/Rickettsia ที่ระดับ genus ได้ยาก จึงใช้ alignment เสริม
- **มี species-aware decontamination filter (V4)** — รวม Burkholderia species-level safeguard, NTC-only removal, ultra-low abundance noise filter และ PC spike-in bypass
- **ผ่าน analytical validation ที่ดี** — 72.1% sample-level analytical performance, 100% specificity, 100% inter-run reproducibility

ผลการศึกษา cohort สรุปได้ว่า

- ~13% ของ cohort มี **Rickettsiales rescue evidence** ซึ่งสอดคล้องกับ endemic epidemiology ของไทย
- ~5% มี ***Mycoplasmopsis*** เป็น candidate signal ของเชื้อ fastidious cell-wall-deficient
- มี ***Leptospira*** 1 ตัวอย่าง และ ***Brucella*** 3 ตัวอย่าง เป็น candidate signals
- **ไม่พบ *B. pseudomallei*** ในระดับ species ใด ๆ (corrected จากเอกสารฉบับก่อน)
- 17% เป็น pre-sequencing failure
- ข้อมูล cohort สนับสนุนว่าเป็น **pilot, hypothesis-generating study** — ทุก candidate ควรยืนยันด้วยวิธีอื่นก่อนรายงานผลทางคลินิก

ทีมห้องปฏิบัติการสามารถนำเนื้อหานี้ไปต่อยอดเขียนต้นฉบับเพิ่มเติมในส่วนที่เป็นข้อมูล wet-lab, clinical, epidemiological และ ethics ที่ทีม bioinformatics ไม่มีข้อมูลโดยตรง ตามรายการที่ section 9 ของเอกสารนี้

---

**หากมีข้อสงสัยเพิ่มเติม หรือต้องการให้อธิบายส่วนใดเพิ่ม ทีม bioinformatics ยินดีให้คำตอบครับ/ค่ะ**

*เอกสารฉบับนี้เขียนวันที่ 2026-05-12 ทุกตัวเลขในเอกสารนี้สามารถยืนยันได้จากไฟล์ `.calls.tsv` และ kreport ของ Centrifuger ใน `/Users/peerahemarajata/Downloads/AFI_P_Final/`*
