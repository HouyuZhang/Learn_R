library(RColorBrewer)

#To display all the color palettes or specific one
display.brewer.all()
display.brewer.all(colorblindFriendly = TRUE)
display.brewer.pal(n = 3, name = "Dark2")

brewer.pal.info
brewer.pal(n = 8, name = "Dark2")


colorRampPalette(brewer.pal(n = 8, name = "Dark2"))(16)

