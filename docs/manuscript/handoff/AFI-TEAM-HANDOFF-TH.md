# เอกสารสรุปการวิเคราะห์ 16S สำหรับทีม AFI (ฉบับย่อ)

**วันที่:** 2026-05-25
**ผู้จัดทำ:** ทีม Bioinformatics (Peera Hemarajata + Yuyi)
**ทีม bioinformatics สนับสนุนได้ถึง:** 2026-05-31

> **วัตถุประสงค์ของเอกสารฉบับนี้**
>
> เอกสารฉบับนี้เป็นเวอร์ชัน**ย่อ** ของเอกสารวิเคราะห์ที่เราเตรียมไว้ เน้นให้ทีมห้องปฏิบัติการ AFI และทีมระบาดวิทยา**อ่านเข้าใจได้รวดเร็ว** โดยไม่ต้องมีพื้นฐาน bioinformatics เอกสารฉบับเต็มภาษาอังกฤษ (`MANUSCRIPT-PACKAGE-FINAL.md`, `MANUSCRIPT-RESULTS-SECTION.md`, `MANUSCRIPT-DISCUSSION-SECTION.md`) อยู่ในโฟลเดอร์ [`docs/manuscript/`](..) (โฟลเดอร์แม่ของไฟล์นี้) ใช้เป็นเอกสารอ้างอิงเมื่อต้องการรายละเอียดเพิ่มเติม
>
> เอกสารคู่กันที่ควรอ่านพร้อมกัน:
> - [`AFI-TEAM-HANDOFF.md`](AFI-TEAM-HANDOFF.md) — ฉบับภาษาอังกฤษ (มีข้อความที่สามารถนำไปใส่ในต้นฉบับโดยตรง)
> - [`LEADS-AND-FRAMING.md`](LEADS-AND-FRAMING.md) — แนวทางการตีกรอบเรื่อง, citation เบื้องต้น, คำถามที่ทีมแพทย์/ระบาดควรตอบ, โครงร่างส่วนที่ทีมต้องเขียนเอง

---

## ส่วนที่ 1 — ปัญหาที่เรากำลังตอบ

โครงการนี้ศึกษาผู้ป่วย Acute Febrile Illness (AFI) คือผู้ป่วยที่มาด้วยอาการ**ไข้เฉียบพลัน** ในภาคตะวันออกเฉียงเหนือของประเทศไทย ผู้ป่วยทุกคนในกลุ่มศึกษาได้รับการเพาะเชื้อจากเลือด (blood culture) ในระหว่างนอนโรงพยาบาล และพบว่า**ขวดเพาะเชื้อให้สัญญาณบวก** จากเครื่องตรวจอัตโนมัติ (เครื่องบอกว่ามีจุลชีพเจริญอยู่ในขวด) แต่เมื่อนำมาเพาะต่อบนวุ้นเลือดในห้องปฏิบัติการ (subculture) **ไม่พบเชื้อขึ้น**ในระยะเวลามาตรฐาน 5–7 วัน โรงพยาบาลที่เก็บตัวอย่างไม่ได้เพาะเชื้อในสภาวะไม่มีออกซิเจน (anaerobic culture) ตามโพรโทคอลของหน่วยงาน

**คำถามวิจัย:** จุลชีพอะไรที่ทำให้ขวดเพาะเชื้อบวก แต่เพาะต่อไม่ขึ้น? การตรวจ 16S rRNA gene สามารถช่วยตอบคำถามนี้ได้ เพราะ 16S ตรวจพบ DNA ของเชื้อได้แม้เชื้อนั้นจะเพาะไม่ขึ้นบนวุ้นปกติ (เช่นเชื้อ Rickettsiales ที่เป็นปรสิตในเซลล์)

**ข้อระวังตั้งแต่ต้น:** ตัวอย่างเลือดที่เราตรวจ 16S **เป็นเลือดผู้ป่วยที่เก็บแยกต่างหาก ไม่ใช่ตัวอย่างจากขวดเพาะเชื้อโดยตรง** ดังนั้นเชื้อที่ตรวจพบใน 16S **อาจ**สัมพันธ์กับสัญญาณบวกในขวด แต่ก็พิสูจน์ความเป็นสาเหตุไม่ได้ ความสัมพันธ์เป็น correlative ไม่ใช่ causal

---

## ส่วนที่ 2 — ขั้นตอนการวิเคราะห์ทั้งหมด 4 ขั้น

เราออกแบบ pipeline การวิเคราะห์ 4 ขั้น เขียนด้วย Workflow Description Language (WDL) รันบน Terra.bio cloud platform โค้ดต้นฉบับเผยแพร่ที่ https://github.com/PHemarajata/afi_terra

### ขั้นที่ 1 — เตรียม reads

นำ paired-end Illumina FASTQ มา (1) ลบลำดับที่เป็น DNA ของมนุษย์ออกด้วยเครื่องมือ NCBI SRA Human Scrubber และ (2) ตัด adapter + กรองคุณภาพด้วยโปรแกรม **fastp** เพื่อให้เหลือเฉพาะ reads ที่สะอาดและไม่ใช่ของคนสำหรับการวิเคราะห์ต่อ

### ขั้นที่ 2 — จำแนกชนิดของแบคทีเรียในระดับสกุล (genus) ด้วยฐานข้อมูลกว้าง

นำ reads ที่สะอาดแล้วมาเทียบกับฐานข้อมูลจีโนมแบคทีเรียและอาร์เคียทั้งหมดของ NCBI โดยเพิ่มเติมจีโนม Rickettsiales เข้าไปด้วย ใช้โปรแกรม **Centrifuger v1.1.0** (เป็น classifier รุ่นใหม่ที่พัฒนาต่อจาก Centrifuge) ผลที่ได้คือจำนวน read ที่ถูกจัดอยู่ในแต่ละสกุล โดยเรารวมตัวเลขในระดับ clade ขึ้นไปถึงระดับ genus

**กฎการประกาศ "Detected"** ในขั้นนี้: เชื้อสกุลใดสกุลหนึ่งจะถูกประกาศว่าตรวจพบก็ต่อเมื่อ **(ก) มีจำนวน read อย่างน้อย 500 และ (ข) จำนวน read ต้องสูงกว่าค่าสูงสุดของ NTC ในรอบเดียวกันอย่างน้อย 5 เท่า** เกณฑ์นี้กรองสัญญาณ background ปกติออกได้พอสมควร

### ขั้นที่ 3 — ค้นหา Rickettsiales อย่างเฉพาะเจาะจง (rescue)

ช่วง V1-V3 ของ 16S **ไม่สามารถแยกระหว่าง Rickettsiales หลายสกุลได้ชัดเจน** (Orientia, Rickettsia, Anaplasma, Ehrlichia) แต่ Rickettsiales เป็นกลุ่มเชื้อ "cannot miss" ที่สำคัญที่สุดในผู้ป่วย AFI ในภาคอีสาน (สครับไทฟัส, สปอตเต็ดฟีเวอร์) เพราะเป็นปรสิตในเซลล์ที่เพาะไม่ขึ้นบนวุ้นปกติ เราจึงเพิ่มขั้น "rescue" โดยนำ reads ไปจัดเรียง (align) กับฐานอ้างอิง Rickettsiales 7 จีโนมที่คัดเลือก ด้วยโปรแกรม **minimap2**

จากนั้นใช้กรอบรายงาน 2 ระดับ:

- **Tier 1 (Confirmed, ระดับสกุล):** mapped reads ≥ 100, breadth ของ coverage ≥ 0.25, และ reads ≥ 5 เท่าของ NTC ในรอบเดียวกัน → รายงานเป็นสกุลเฉพาะ (Orientia หรือ Rickettsia)
- **Tier 2 (Probable, ระดับ order):** mapped reads ≥ 50, breadth ≥ 0.20, และ reads สูงกว่า NTC แต่ไม่ถึง 5 เท่า → รายงานเป็น **"Rickettsiales detected (ไม่สามารถระบุสกุลแน่ชัด, แนะนำตรวจ qPCR ยืนยัน)"** ในบริบทคลินิกของภาคอีสาน Tier 2 ก็เพียงพอที่จะพิจารณาให้ยา doxycycline ตามอาการอยู่แล้ว

### ขั้นที่ 4 — กรอง contamination (V4 filter)

เลือดเป็นตัวอย่างที่มีแบคทีเรียน้อยมาก (low-biomass specimen) ดังนั้นสัญญาณ false-positive จากการปนเปื้อนของน้ำยาและชุดทดสอบเป็นปัญหาที่ทราบกันดี (Salter et al., 2014) เราสร้าง filter 4 ชั้น (เรียกว่า V4 filter) เพื่อตัดสัญญาณปนเปื้อนออกโดยไม่ให้ลบเชื้อก่อโรคที่แท้จริงไปด้วย

**ประเด็นที่ทีมต้องเข้าใจ 3 ข้อ:**

1. **กรองเชื้อ kit/skin/water contaminant อย่างเข้มงวด** — มีรายการ 11 สกุลที่ถูกลบทันทีเมื่อพบ เช่น *Cutibacterium*, *Staphylococcus*, *Brevundimonas*, *Acinetobacter*, *Corynebacterium*, *Pseudomonas* (รายการเต็มอยู่ใน [`docs/manuscript/scripts/afi_decontamination_filter_v4.py`](../scripts/afi_decontamination_filter_v4.py))
2. **Burkholderia มี species-level safeguard** — สกุล Burkholderia มีทั้ง *B. pseudomallei* (เมลิออยโดสิส) และ *B. cepacia complex* (ปนเปื้อนจากน้ำยา) ดังนั้นเมื่อเจอสัญญาณ Burkholderia ระบบจะอ่าน Centrifuger report ในระดับ species ใหม่: จะเก็บไว้เป็น *B. pseudomallei* เฉพาะเมื่อ (ก) มี species reads ≥ 500 และ (ข) จำนวน species reads สูงกว่า NTC ในรอบเดียวกัน
3. **PC spike-in bypass** — สำหรับตัวอย่าง positive control การลบ Tier-A จะถูก skip สำหรับเชื้อที่เป็น spike-in ที่ควรพบในตัวอย่างนั้น (เช่น *Pseudomonas* ใน PC_SINGLE *P. aeruginosa*)

---

## ส่วนที่ 3 — ผลการตรวจสอบความถูกต้องของ pipeline (validation panel)

ใช้ตัวอย่าง 48 ตัวอย่าง (33 คลินิก + 10 positive control + 5 NTC) ครอบคลุม 9 รอบ sequencing

| ประเภท | ถูกต้อง / รวม | อัตรา |
|---|---|---|
| Clinical (เข้มงวด) | 21 / 33 | **63.6%** |
| Positive control | 10 / 10 | **100%** |
| Negative template control (specificity) | 5 / 5 | **100%** |
| **Sample-level analytical (clinical + PC)** | **31 / 43** | **72.1%** |
| Overall validation accuracy | 36 / 48 | 75.0% |

**Reproducibility ระหว่างรอบ (inter-run): 100%** ทุก control ผ่านเกณฑ์ในทุกรอบ

ตัวอย่าง clinical ที่ไม่ตรงกับเป้า 12 ตัวอย่างมี 2 รูปแบบ:
- **2 ตัวอย่าง failed pre-sequencing** (อ่านไม่ได้เลย, 0 reads) เป็นปัญหาขั้น DNA extraction / library prep / sequencing depth ไม่ใช่ pipeline error
- **10 ตัวอย่าง abundance-driven discordance** คือ pipeline เห็นเชื้อเป้าหมายอยู่ แต่มีเชื้ออื่นเด่นกว่าในตัวอย่างนั้น (สะท้อนธรรมชาติของ 16S ในตัวอย่างที่มีจุลชีพหลายชนิด)

ตัวอย่าง validation panel ของ *B. pseudomallei* 3 ตัวอย่าง (`09502813_S2_L001`, `09-0-02165`, `09700912_S3_L001`) มี species-level reads 65,016 / 20,061 / 13,744 ตามลำดับ ผ่าน safeguard ได้สบาย → ยืนยันว่า pipeline ตรวจเมลิออยโดสิสได้เมื่อมีปริมาณเชื้อจริง

---

## ส่วนที่ 4 — ผลในกลุ่มศึกษา (86 ตัวอย่าง)

ตัวอย่าง 86 ตัวอย่างกระจายอยู่ใน 5 รอบ sequencing ทั้งหมดมีไฟล์ output

- **71 ตัวอย่าง (82.6%) มี ≥1 positive call** (Detected / Confirmed / Probable)
- **15 ตัวอย่าง (17.4%) ไม่พบเชื้อเลย** (เข้าได้กับ pre-sequencing failure)

V4 filter ถูกใช้กับ 217 genus-detection ทั้งหมดในกลุ่มที่มี signal → **ลบออก 70 (32.3%)** เก็บไว้ 147

### 4.1 Rickettsiales rescue — ผลลัพธ์สำคัญที่สุด

**พบหลักฐาน Rickettsiales rescue ใน 11 ของ 86 ตัวอย่าง (12.8%; เทียบเท่า 11 ของ 71 ตัวอย่างที่มีการตรวจพบ = 15.5%)** ประกอบด้วย **1 ตัวอย่าง Tier-1 genus-level Orientia call** (ตัวอย่าง `16901195_S5_L001`, breadth 0.3235, NTC headroom ~11.5 เท่า — เป็น Rickettsiales call ที่ชัดที่สุดในกลุ่ม) และ **10 ตัวอย่าง Tier-2 order-level rescue** (breadth 0.21–0.24, อยู่ใต้เกณฑ์ Tier-1 เล็กน้อย)

ทั้ง 11 ตัวอย่างควรส่งทำ Rickettsiales-specific qPCR เพื่อยืนยันชนิด; การให้ doxycycline เป็น empirical ในบริบทนี้เป็นแนวปฏิบัติมาตรฐานอยู่แล้ว ตารางรายตัวอย่างพร้อม breadth + NTC + confidence ดูได้ใน [`AFI-TEAM-HANDOFF.md`](AFI-TEAM-HANDOFF.md) §D.1 หรือ [`docs/manuscript/appendices/APPENDIX-STUDY-SAMPLES.md`](../appendices/APPENDIX-STUDY-SAMPLES.md)

### 4.2 ไม่พบเมลิออยโดสิสในกลุ่มศึกษา

มีสัญญาณ Burkholderia ที่ระดับ genus ใน 5 ตัวอย่าง แต่เมื่ออ่าน Centrifuger ระดับ species ใหม่ **ไม่มีตัวอย่างใดมี *B. pseudomallei* reads ถึงเกณฑ์ 500 และไม่มีตัวอย่างใดที่ species reads สูงกว่า NTC** สัญญาณ genus ในตัวอย่างเหล่านี้มาจาก *B. cepacia complex* ซึ่งเป็นเชื้อปนเปื้อนจาก kit/น้ำ V4 filter ลบสัญญาณเหล่านี้ออกอย่างถูกต้องทั้ง 5 ตัว

**ข้อแก้ไขจากฉบับก่อนหน้า:** รายงานก่อนหน้านี้ระบุว่า `23200430_S6_L001` มี *B. pseudomallei* 54,126 reads ที่ 15.56% นั่นเป็นข้อผิดพลาด — จำนวนนั้นเป็น **genus-level count** ไม่ใช่ species-level จำนวน *B. pseudomallei* ที่ระดับ species ในตัวอย่างนี้คือ 8 reads ซึ่งต่ำกว่า NTC ของรอบเดียวกัน (10–51 reads)

### 4.3 Mycoplasmopsis — สัญญาณกลุ่ม fastidious ที่เด่นที่สุด

พบ *Mycoplasmopsis* ใน **4 ของ 86 ตัวอย่าง (4.7%)** mean abundance 39.78%, max 70.55% เชื้อกลุ่ม Mycoplasma ไม่มีผนังเซลล์ ต้องใช้อาหารเลี้ยงเชื้อพิเศษ (เติม sterol) และโตช้า (1–3 สัปดาห์) จึงไม่น่าจะขึ้นบนวุ้นเลือดธรรมดาในเวลามาตรฐาน 5–7 วัน → สัญญาณนี้สอดคล้องกับ phenotype "ขวดบวก แต่เพาะต่อไม่ขึ้น" แต่จำนวนตัวอย่างน้อย ยังไม่สามารถสรุปได้ในระดับ cohort → ควรทำ Mycoplasma-specific PCR ยืนยันก่อนใช้ทางคลินิก

### 4.4 Leptospira (1 ตัวอย่าง) และ Brucella (3 ตัวอย่าง)

ตัวอย่าง `09801652_S5_L001` มี *Leptospira* 7,294 reads ที่ 22.78% **แต่** ในรอบ sequencing เดียวกัน NTC (`NTC2_ExDw_S13_L001`) มี *Leptospira* 78,691 reads (สูงกว่าตัวอย่างศึกษาประมาณ 10 เท่า) pipeline ดูเหมือนจะใช้ NTC ที่ไม่ปนเปื้อนในการคำนวณ NCmax สำหรับตัวอย่างนี้ แต่ยังคงเป็น flag ที่ต้องระวัง → ควรยืนยันด้วย serology + Leptospira-specific qPCR ก่อนสรุปทางคลินิก

*Brucella* ใน 3 ตัวอย่างที่ความเข้มข้นต่ำมาก (mean 0.98%, max 1.96%) อยู่ใกล้กับ noise floor → ต้องยืนยันด้วย serology + Brucella-specific qPCR

---

## ส่วนที่ 5 — ข้อจำกัดสำคัญที่ห้ามลืม

ข้อจำกัด 4 ข้อนี้**ต้องคงไว้ในต้นฉบับ** อย่าตัดออกหรือลดความหนักแน่นลง

1. **Sample-source mismatch** — ตัวอย่าง 16S เป็นเลือดผู้ป่วย ไม่ใช่ตัวอย่างจากขวดเพาะเชื้อ การตรวจพบเชื้อใน 16S สอดคล้องกับสาเหตุของสัญญาณบวกในขวด แต่พิสูจน์ความเป็นสาเหตุไม่ได้
2. **ไม่ประเมิน viability** — 16S ตรวจ DNA จากทั้งเซลล์มีชีวิต เซลล์ตาย และเซลล์ที่ viable-but-non-culturable (VBNC) ดังนั้นสัญญาณสูงไม่ได้แปลว่าเชื้อยังมีชีวิตหรือเพาะขึ้นได้
3. **V1-V3 primer biases** — primer 27F มีอคติในการตรวจจับ Gram-positive anaerobes บางกลุ่ม, Mycobacterium, *Bifidobacterium*/*Gardnerella*/*Atopobium* บางสปีชีส์ → "ไม่พบในการตรวจนี้" ≠ "ไม่มีในตัวอย่างจริง"
4. **NTC contamination ในรอบ 6_and_7** — มี NTC ที่ปนเปื้อน *Leptospira*, *Burkholderia*, *Brevundimonas* จำนวนมาก ตรงนี้ต้องระบุไว้เป็น caveat ของ Leptospira candidate ใน `09801652_S5_L001`

ข้อจำกัดเพิ่มเติม 6 ข้อ (ขนาด pilot, ข้อจำกัดความละเอียดระดับ species, การ depend on dataset ของ filter, ไม่มี clinical outcome ในแพ็กเกจ bioinformatics, ไม่มี anaerobic culture เปรียบเทียบ, และ practical detection floor ของ *B. pseudomallei*) อยู่ใน back-pocket §5 ([`MANUSCRIPT-PACKAGE-FINAL.md`](../MANUSCRIPT-PACKAGE-FINAL.md))

---

## ส่วนที่ 6 — ใครรับผิดชอบอะไร และจะหาเอกสารฉบับเต็มได้ที่ไหน

**ทีม Bioinformatics รับผิดชอบ:**
- §2 Methods (analysis pipeline + concordance definitions) → ใช้ §B ใน [`AFI-TEAM-HANDOFF.md`](AFI-TEAM-HANDOFF.md)
- §3 Results (validation + study cohort) → ใช้ §C + §D ใน [`AFI-TEAM-HANDOFF.md`](AFI-TEAM-HANDOFF.md)
- §4 Discussion (4 ย่อหน้าตีความข้อมูลของเรา) → ใช้ §E ใน [`AFI-TEAM-HANDOFF.md`](AFI-TEAM-HANDOFF.md)
- ส่วนที่เกี่ยวกับ Methods + Results ใน Abstract → ใช้ §F ใน [`AFI-TEAM-HANDOFF.md`](AFI-TEAM-HANDOFF.md)
- ตารางรายตัวอย่าง (appendix) → อยู่ใน [`docs/manuscript/appendices/`](../appendices/)
- รูปที่ 1 (Sankey) → อยู่ใน [`docs/manuscript/figures/`](../figures/)
- Reference สำหรับซอฟต์แวร์ → อยู่ใน [`docs/manuscript/MANUSCRIPT-PACKAGE-FINAL.md`](../MANUSCRIPT-PACKAGE-FINAL.md) §8

**ทีมห้องปฏิบัติการ (NIH wet-lab) รับผิดชอบ:**
- ส่วน Background / Introduction ของ Abstract
- §2 Methods: sample collection, DNA extraction kit, V1-V3 primer sequences, library prep kit, MiSeq configuration
- ส่วน Discussion ที่เกี่ยวกับ wet-lab methodology
- รายชื่อผู้แต่ง + funding + CoI (ประสานกับ CDC)

**ทีมระบาด (EPI) รับผิดชอบ:**
- §1 Introduction / Background (epidemiology ของ AFI ในภาคอีสาน)
- §2 Methods: IRB approval, enrolment period, hospital/district, inclusion/exclusion criteria
- §4 Discussion ส่วน clinical context, case vignettes (per-case treatment, response, exposure, outcome), clinical implications
- §6 Conclusion
- Reference สำหรับ AFI epidemiology, Rickettsioses, melioidosis, leptospirosis, brucellosis → ดู starter list ใน [`LEADS-AND-FRAMING.md`](LEADS-AND-FRAMING.md)

**ทีม Bioinformatics สนับสนุนได้ถึง 2026-05-31** หลังจากนั้นจะไม่สามารถตอบคำถามได้ทันที กรุณาอ่านเอกสารฉบับนี้ + เอกสารฉบับเต็มก่อนถามคำถาม ถ้ามีคำถามที่ตอบไม่ได้จากเอกสาร ส่งภายในวันที่ 30 พฤษภาคมเพื่อให้มีเวลาตอบ
