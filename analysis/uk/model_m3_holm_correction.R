m1 <- read.csv("results/uk/primary/M1_interaction_coefficient.csv")
m2 <- read.csv("results/uk/primary/M2_interaction_coefficient.csv")

raw_p <- c(m1$p_value[1], m2$p_value[1])
adj_p <- p.adjust(raw_p, method="holm")

out <- data.frame(
  model=c("M1_manage1","M2_manage2"),
  interaction=c(as.character(m1$term[1]), as.character(m2$term[1])),
  estimate=c(m1$estimate[1], m2$estimate[1]),
  raw_p=raw_p,
  holm_p=adj_p,
  reject_0_05=adj_p < 0.05
)

write.csv(out, "results/uk/primary/M3_holm_moderation_family.csv", row.names=FALSE)

cat("M3_HOLM_PASS\n")
print(out)
