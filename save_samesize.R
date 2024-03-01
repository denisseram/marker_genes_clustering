library(dplyr)
library(ggplot2)

# CONVERTIR ARI
df_ari <- data.frame(ARI_matrix)

new_ari <- df_ari[2]
new_ari$dataset <- "baron.mouse_markergenes"
new_ari$algorithm <- rownames(new_ari)
rownames(new_ari) <- NULL
names(new_ari)[1] <- "ARI"

new_ari1 <- df_ari[1]
new_ari1$dataset <- "baron.mouse"
new_ari1$algorithm <- rownames(new_ari1)
rownames(new_ari1) <- NULL
names(new_ari1)[1] <- "ARI"

merged_df <- rbind(new_ari, new_ari1)
merged_df$ARI <- as.numeric(as.character(merged_df$ARI))

ari_plot <- ggplot(data = merged_df, aes(x = algorithm, y = ARI, fill = dataset)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_manual(values = c("#EABC40", "#72B1BF"))

ggsave("ari_datasets_algorithms.png", plot = ari_plot, width = 8, height = 6)

# CONVERTIR AMI
df_ami <- data.frame(AMI_matrix)

new_ami <- df_ami[2]
new_ami$dataset <- "baron.mouse_markergenes"
new_ami$algorithm <- rownames(new_ami)
rownames(new_ami) <- NULL
names(new_ami)[1] <- "AMI"

new_ami1 <- df_ami[1]
new_ami1$dataset <- "baron.mouse"
new_ami1$algorithm <- rownames(new_ami1)
rownames(new_ami1) <- NULL
names(new_ami1)[1] <- "AMI"

merged_df <- rbind(new_ami, new_ami1)
merged_df$AMI <- as.numeric(as.character(merged_df$AMI))

ami_plot <- ggplot(data = merged_df, aes(x = algorithm, y = AMI, fill = dataset)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_manual(values = c("#EABC40", "#72B1BF"))

ggsave("ami_datasets_algorithms.png", plot = ami_plot, width = 8, height = 6)

# CONVERTIR VI
df_vi <- data.frame(VI_matrix)

new_vi <- df_vi[2]
new_vi$dataset <- "baron.mouse_markergenes"
new_vi$algorithm <- rownames(new_vi)
rownames(new_vi) <- NULL
names(new_vi)[1] <- "VI"

new_vi1 <- df_vi[1]
new_vi1$dataset <- "baron.mouse"
new_vi1$algorithm <- rownames(new_vi1)
rownames(new_vi1) <- NULL
names(new_vi1)[1] <- "VI"

merged_df <- rbind(new_vi, new_vi1)
merged_df$VI <- as.numeric(as.character(merged_df$VI))

vi_plot <- ggplot(data = merged_df, aes(x = algorithm, y = VI, fill = dataset)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_y_continuous(limits = c(0, 5)) +
  scale_fill_manual(values = c("#EABC40", "#72B1BF"))

ggsave("vi_datasets_algorithms.png", plot = vi_plot, width = 8, height = 6)


# CONVERTIR ARI
df_ari <- data.frame(number_clusters)

new_ari <- df_ari[2]
new_ari$dataset <- "baron.mouse_markergenes"
new_ari$algorithm <- rownames(new_ari)
rownames(new_ari) <- NULL
names(new_ari)[1] <- "Number_of_clusters"

new_ari1 <- df_ari[1]
new_ari1$dataset <- "baron.mouse"
new_ari1$algorithm <- rownames(new_ari1)
rownames(new_ari1) <- NULL
names(new_ari1)[1] <- "Number_of_clusters"

merged_df <- rbind(new_ari, new_ari1)
merged_df$Number_of_clusters <- as.numeric(as.character(merged_df$Number_of_clusters))

num_p <- ggplot(data = merged_df, aes(x = algorithm, y = Number_of_clusters, fill = dataset)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_y_continuous(limits = c(0, 206)) +
  scale_fill_manual(values = c("#EABC40", "#72B1BF"))

num_p <- ggplot(data = merged_df, aes(x = algorithm, y = Number_of_clusters, fill = dataset)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_hline(yintercept = 13, linetype = "dashed", color = "red") +  # Agregar línea horizontal en y = 13
  scale_y_continuous(limits = c(0, 206)) +
  scale_fill_manual(values = c("#EABC40", "#72B1BF"))


ggsave("clu_num_datasets_algorithms.png", plot = num_p, width = 8, height = 6)
