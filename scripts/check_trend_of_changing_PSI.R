#Use the data from "all_AS_events_TPM.R" to check whether PI_H and PI_L change in different ways in sig-AS events or non-sig AS events

all_AS_events_bed_TPM_q05_raw %>%
  filter(comp == "HS1.vs.HS6") %>%
  filter(PI == "PI_H" | PI == "PI_L") %>%
  mutate(Inc_1H = (Inclvl1_r1 + Inclvl1_r2)/2, Inc_6H = (Inclvl2_r1 + Inclvl2_r2)/2) %>%
  mutate(Inc_rise = if_else(Incdf > 0, "yes", if_else(Incdf < 0, "no", "equal"))) %>%
  group_by(AS_type, PI, Inc_rise) %>%
  summarise(count = n()) %>%
  View()


test <- all_DAS_events_bed_TPM_q05_raw %>%
  filter(comp == "HS1.vs.HS6") %>%
  filter(PI == "PI_H" | PI == "PI_L") %>%
  mutate(Inc_1H = (Inclvl1_r1 + Inclvl1_r2)/2, Inc_6H = (Inclvl2_r1 + Inclvl2_r2)/2) %>%
  ggplot() +
  geom_point(aes(x = Inc_1H, y = Inc_6H)) +
  facet_wrap(~ AS_type)


all_DAS_events_bed_TPM_q05_raw %>%
  filter(comp == "HS1.vs.HS6") %>%
  filter(PI == "PI_H" | PI == "PI_L") %>%
  mutate(Inc_1H = (Inclvl1_r1 + Inclvl1_r2)/2, Inc_6H = (Inclvl2_r1 + Inclvl2_r2)/2) %>%
  mutate(Inc_rise = if_else(Incdf > 0, "yes", if_else(Incdf < 0, "no", "equal"))) %>%
  group_by(comp, AS_type, PI, Inc_rise) %>%
  summarise(count = n()) %>%
  View()

all_DAS_events_bed_TPM_q05_raw %>%
  filter(comp == "HS6.vs.HS0") %>%
  filter(PI == "PI_H" | PI == "PI_L") %>%
  mutate(Inc_6H = (Inclvl1_r1 + Inclvl1_r2)/2, Inc_0H = (Inclvl2_r1 + Inclvl2_r2)/2) %>%
  mutate(Inc_rise = if_else(Incdf < 0, "yes", if_else(Incdf > 0, "no", "equal"))) %>%
  group_by(comp, AS_type, PI, Inc_rise) %>%
  summarise(count = n()) %>%
  View()
