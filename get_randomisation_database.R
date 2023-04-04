## Load randomisation data to get ITT database for the molnupiravir final analysis


data.TH1 <- read.csv("~/Dropbox/PLATCOV/data-TH1.csv")
data.TH1$Date = as.POSIXct(data.TH1$Date,format='%a %b %d %H:%M:%S %Y')
data.TH1 = data.TH1[data.TH1$Date > '2022-06-05' & data.TH1$Date <= '2023-02-20',  ]

data.TH1$ID = paste('PLT-TH1-',data.TH1$randomizationID,sep='')

xx = data.TH1[, c('ID', 'Treatment')]

library(stringr)
for(i in 1:nrow(xx)){
  id = unlist(strsplit(xx$ID[i],split = '-'))
  id[3] = str_pad(id[3], 3, pad = "0")
  id = paste(id, collapse = '-')
  xx$ID[i]=id
}

table(xx$Treatment)

write.csv(x = xx, file = 'ITT_population.csv')

