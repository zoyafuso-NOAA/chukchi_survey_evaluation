load(paste0("presentations/2022/WKUSER_2022/results/Str_3/",
            "random_survey_sim_results.RData"))
load(paste0("presentations/2022/WKUSER_2022/results/Str_3/",
            "sys_survey_sim_results.RData"))

# load(paste0("results/chukchi_otter/sys_survey_sim_results.RData"))


png(filename = "presentations/2022/WKUSER_2022/results/True_CV.png", 
    width = 8, height = 6, units = "in", res = 200)
par(mfrow = c(2, 3), mar = c(3, 5, 1, 1), oma = c(3, 3, 0, 1))
spp_list <- dimnames(true_cv_ms_strs)[[2]]

for (ispp in spp_list) {
  plot_this <- cbind(true_cv_srs[,ispp], true_cv_ms_strs[, ispp], true_cv_sys[, ispp])
  matplot(x = target_n, y = plot_this,
          las = 1, type = "l", lty = 1, lwd = 2, cex.axis = 1.5,
          ylab = "", xlab = "", ylim = c(0, max(plot_this)),
          col = c("black", "red", "blue"))
  legend("topright", legend = ispp, bty = "n", cex = 1.5, text.font = 2)
}
mtext(side = 1, text = "Sample Size", line = 1, font = 2, outer = TRUE, cex = 2)
mtext(side = 2, text = "True CV", line = 0.5, font = 2, outer = TRUE, cex = 2)
dev.off()

png(filename = "presentations/2022/WKUSER_2022/results/RRMSE_CV.png", 
    width = 8, height = 6, units = "in", res = 200)
par(mfrow = c(2, 3), mar = c(3, 5, 1, 1), oma = c(3, 3, 0, 1))
spp_list <- dimnames(true_cv_ms_strs)[[2]]
for (ispp in spp_list) {
  plot_this <- cbind(rrmse_cv_srs[,ispp], rrmse_cv_ms_strs[, ispp], rrmse_cv_sys[, ispp])
  matplot(x = target_n, y = plot_this,
          las = 1, type = "l", lty = 1, lwd = 2, cex.axis = 1.5,
          ylab = "", xlab = "", ylim = c(0, 1.1*max(plot_this)),
          col = c("black", "red", "blue"))
  legend("topright", legend = ispp, bty = "n", cex = 1.5, text.font = 2)
}
mtext(side = 1, text = "Sample Size", line = 1, font = 2, outer = TRUE, cex = 2)
mtext(side = 2, text = "RRMSE", line = 0.5, font = 2, outer = TRUE, cex = 2)
dev.off()

png(filename = "presentations/2022/WKUSER_2022/results/legend.png", 
    width = 10, height = 1, units = "in", res = 200)
par(mfrow = c(1, 1), mar = c(0, 0, 0, 0), oma = c(0,0,0,0))
plot(1, type = "n", axes = F)
legend("center", lty = 1, lwd = 3, col = c("black", "red", "blue"), 
       cex = 1, bty = "n", ncol = 3, 
       legend = c("Simple\nRandom", "Stratified", "Systematic"))
dev.off()
