# given a directory containing networks for each sample, this script reads the data

networks.dir = "./metabolic_nets/"

cat(paste("reading network data from ", networks.dir, "\n", sep=""))

network.files = list.files(path=networks.dir, pattern="*.net")

sample.names = sapply( strsplit(as.character(network.files), '\\.'), "[", 1) # note '\\.' to match . character and sapply and "[" for pulling the first column
num.samples = length(sample.names)
for (i in 1:num.samples) { sample.names[i] = substr(network.files[i], 1, nchar(network.files[i]) - 14) }

excludedSamples = c() 


S = list();
R = list();
Reactions = list();
Compounds = c()
cur.compound = 0

Compounds.names = c();
cat("collecting compound sets\n")
for (network.file in network.files){
  if (network.file %in% excludedSamples){
    cat(paste( "skipping excluded",  network.file ,"\n", sep=" "))
    next
  }
  cat(paste( "collecting compounds in ",  network.file , "\n", sep=" "))
  lines = read.table( paste(networks.dir, network.file, sep=""), sep=",", blank.lines.skip=FALSE, fill=TRUE)
  Compounds.names = c(Compounds.names, as.character(lines[,1]), as.character(lines[,3]) )
}
#browser()
Compounds.names = unique(Compounds.names);
Compounds.names = setdiff( Compounds.names, c("") )
Compounds = seq(0, length(Compounds.names)-1 )
names(Compounds) = Compounds.names
cat(paste( length(Compounds) , " compounds in all networks \n", sep=""))

cat("reading networks \n")
for (network.file in network.files){
  if (network.file %in% excludedSamples){
    cat(paste( "skipping excluded",  network.file ,"\n", sep=" "))
    next
  }
  cat(paste( "processing",  network.file ,"\n", sep=" "))
  lines = read.table( paste(networks.dir, network.file, sep=""), sep=",", blank.lines.skip=FALSE, fill=TRUE)
  lines.count = nrow(lines)
  S[[network.file]] = list()
  R[[network.file]] = list()
  Reactions[[network.file]] = list()

  reaction.borders = which(lines[,1] == "")
  reaction.right.borders = reaction.borders -1
  reaction.left.borders = c(1, reaction.borders[1:(length(reaction.borders)-1)]+1)

  #browser()
  cat(paste("adding ", length(reaction.left.borders), " reactions to network for ", network.file, "\n", sep=""))
  length(S[[network.file]]) = length(reaction.left.borders)
  for ( i in 1:length(reaction.left.borders)){
    S[[network.file]][[i]] = Compounds[ as.vector( as.character(lines[(reaction.left.borders[i]):(reaction.right.borders[i]), 1]) ) ]
    Reactions[[network.file]][[i]] = unique( as.vector( as.character(lines[(reaction.left.borders[i]):(reaction.right.borders[i]), 2]) ) )
    R[[network.file]][[i]] = Compounds[ as.vector( as.character(lines[(reaction.left.borders[i]):(reaction.right.borders[i]), 3]) ) ]
    
  }
  

}

save.image(file="./results/network.RData")
