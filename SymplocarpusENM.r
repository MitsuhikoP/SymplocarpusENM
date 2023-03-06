library(tidyverse)
library(dismo)
library(rJava)
library(rgeos)
library(usdm)

#samples
samples=read.csv("SapleLocation.csv")
smSre=filter(samples, SP=="Sre")
smSni=filter(samples, SP=="Sni")
colnames(smSre)[2:3]=c("Longitude","Latitude")
colnames(smSni)[2:3]=c("Longitude","Latitude")


###environmental variables
env.files=list.files(path="./ENV_DIR", full.names = TRUE)
env1p4=stack(env.files)
crs(env1p4) <- "+proj=longlat +datum=WGS84 +no_defs"
names(env1p4) <- c("bio01", "bio02", "bio03", "bio04", "bio05", "bio06", "bio07", "bio08", "bio09", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17", "bio18", "bio19")
env1p4 <- setMinMax(env1p4)

env.files=list.files(path="./ENV_DIR_CCMH", full.names = TRUE)
envCCMH=stack(env.files)
crs(envCCMH) <- "+proj=longlat +datum=WGS84 +no_defs"
names(envCCMH) <- c("bio01", "bio02", "bio03", "bio04", "bio05", "bio06", "bio07", "bio08", "bio09", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17", "bio18", "bio19")
envCCMH <- setMinMax(envCCMH)

env.files=list.files(path="./ENV_DIR_MRLGM", full.names = TRUE)
envMRLGM=stack(env.files)
crs(envMRLGM) <- "+proj=longlat +datum=WGS84 +no_defs"
names(envMRLGM) <- c("bio01", "bio02", "bio03", "bio04", "bio05", "bio06", "bio07", "bio08", "bio09", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17", "bio18", "bio19")
envMRLGM <- setMinMax(envMRLGM)

###VIF around sampling points###
cr=circles(samples[,2:3],d=500000,lonlat=TRUE)
pol=gUnaryUnion(cr@polygons)
maskc=mask(env1p4$bio01,pol)
SampleRandom=sampleRandom(x=maskc, size=10000,na.rm=TRUE, ext=pol, xy=TRUE)[,1:2]
randRaster=data.frame(raster::extract(env1p4,rbind(as.matrix(samples[,2:3]), as.matrix(SampleRandom))))
envif=vifstep(randRaster,th=5, maxobservations=11000)
env2=env1p4[[envif@results$Variables]]


###Maxent###############
absPoints=sampleRandom(x=maskc, size=10000, na.rm=TRUE, ext=pol, xy=TRUE)[,1:2]
colnames(absPoints)=c("Longitude","Latitude")
env2Sre = data.frame(cbind(pa=c(rep(1, nrow(smSre)), rep(0, nrow(absPoints))), raster::extract(env2, rbind(smSre[,2:3], absPoints))))
env2Sni = data.frame(cbind(pa=c(rep(1, nrow(smSni)), rep(0, nrow(absPoints))), raster::extract(env2, rbind(smSni[,2:3], absPoints))))

mxSre=maxent(env2Sre[,-1], env2Sre[,1])
mxSni=maxent(env2Sni[,-1], env2Sni[,1])

e_mxSre=evaluate(env2Sre[env2Sre$pa==1,-1], env2Sre[env2Sre$pa==0,-1], mxSre)
e_mxSni=evaluate(env2Sni[env2Sni$pa==1,-1], env2Sni[env2Sni$pa==0,-1], mxSni)

pdf("response_curve.pdf", width=10,height=6)
response(mxSre)
response(mxSni)
dev.off()


predSreCurrent=predict(mxSre, env2)
predSniCurrent=predict(mxSni, env2)
predSreCCMH=predict(mxSre, envCCMH)
predSniCCMH=predict(mxSni, envCCMH)
predSreMRLGM=predict(mxSre, envMRLGM)
predSniMRLGM=predict(mxSni, envMRLGM)


###visualize##
set_map=function(RASTER, SP, ENV){
	mp=data.frame(rasterToPoints(RASTER))
	mp$species=SP
	mp$env=ENV
	colnames(mp)=c("longitude","latitude","Prob","Species","Env")
	return(mp)
}

mpFIG=data.frame()
mpFIG=rbind(mpFIG, set_map(mean(predSreCurrent), "Sre", "Current"))
mpFIG=rbind(mpFIG, set_map(mean(predSniCurrent), "Sni", "Current"))
mpFIG=rbind(mpFIG, set_map(mean(predSreCCMH), "Sre", "MH"))
mpFIG=rbind(mpFIG, set_map(mean(predSniCCMH), "Sni", "MH"))
mpFIG=rbind(mpFIG, set_map(mean(predSreMIROC), "Sre", "LGM"))
mpFIG=rbind(mpFIG, set_map(mean(predSniMIROC), "Sni", "LGM"))

mpFIG$Env = factor( mpFIG$Env, levels = c( "Current", "MH", "LGM" ) ) 

g = ggplot(mpFIG,aes(x=longitude,y=latitude, fill=Prob))
g = g+theme_bw() + theme(panel.background = element_rect(fill = "transparent",color = NA),plot.background = element_rect(fill = "transparent",color = NA) ) 
g = g+coord_fixed(ratio=110977/89011) #latitude 37 degrees north.
g = g+geom_tile()+facet_grid(Species~Env)+xlab("Longitude")+ylab("Latitude")
gv = g+scale_fill_viridis_c()
ggsave(gv, file="plot_map.pdf", width=12,height=9)



###ENMTools
library(ENMTools)

raster.breadth(predSreCurrent)$B1
raster.breadth(predSniCurrent)$B1
raster.breadth(predSreCCMH)$B1
raster.breadth(predSniCCMH)$B1
raster.breadth(predSreMRLGM)$B1
raster.breadth(predSniMRLGM)$B1

raster.overlap(predSreCurrent, predSreCCMH)$D
raster.overlap(predSreCCMH, predSreMRLGM)$D
raster.overlap(predSreCurrent, predSreCCMH)$I
raster.overlap(predSreCCMH, predSreMRLGM)$I
raster.overlap(predSniCurrent, predSniCCMH)$D
raster.overlap(predSniCCMH, predSniMRLGM)$D
raster.overlap(predSniCurrent, predSniCCMH)$I
raster.overlap(predSniCCMH, predSniMRLGM)$I

env3=check.env(env2)
enmSre=enmtools.species(species.name="Srenifolius", presence.points=smSre, background.points=data.frame(SampleRandom))
enmSni=enmtools.species(species.name="Snipponicus", presence.points=smSni, background.points=data.frame(SampleRandom))

id.mx = identity.test(species.1=enmSre, species.2=enmSni,env=env3, type="mx", nreps=999)
bg.mx.sym = background.test(species.1=enmSni, species.2=enmSre, env=env3, type="mx", nreps=999, test.type="symmetric")

