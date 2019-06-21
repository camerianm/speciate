#!/usr/bin/env Rscript

#Download libraries
library(CHNOSZ)
library(dplyr)
library(stringr)
library(tools)

#Get databases to cross-reference
media <- read.csv("from_DG_media_molar_final.csv", header = T)
cpds <- read.delim("ModelSEED_compounds.tsv", header = T)
conversion <- read.csv("CHNOSZ_SEED_key.csv", header = T)


### ID compounds actually present in media
##### (this part takes the longest, maybe around 30 seconds to 2 minutes)
print("Linking cpd ID numbers to compound names. One moment...")
cpds_used <- data.frame()
for (cpdID in names(media)) {
    if (cpdID %in% cpds$id) {
    i <- filter(cpds, str_detect(id, cpdID)==TRUE)
    cpds_used <- rbind(i, cpds_used)
    }
}
CHNOSZ_conv <- merge(cpds_used, conversion, by.x="abbreviation", by.y="SEED_abbrev")
print("done!")

### Identify which columns we want to avoid doing calculations for (so we can exclude any media that contain them).
###### Also, print out filtering commands. (Hacky, sure... but works well enough.)
print("Checking which compounds are in CHNOSZ...")
not_used <-  list()
media_usable <- media
for (cpd in names(media)){
    if (all(str_detect(CHNOSZ_conv$id, cpd)==FALSE) && (str_detect(cpd,'cpd'))){
        not_used <- c(not_used, cpd)
    }
}
filter_command <- paste('media_used <- filter(media_usable, ', paste(not_used, collapse='==0 , '), '==0)')
cols_to_remove <- paste('filtered_media <- select(filtered_media, -c(', paste(not_used, collapse=', '), '))')
#filter_command
#cols_to_remove
print("done!")

### Reduce our media list to only the media and (preliminary/first-pass/query-centric) species we need.
##### ...by copying and pasting what prints from filter_command and cols_to_remove Again, hacky, but works well enough.
#First, we remove from consideration the media (rows) that require compounds missing from CHNOSZ.
print('Removing media which require compounds missing from CHNOSZ...')
filtered_media <- filter(media_usable,  cpd00010==0 , cpd00015==0 , cpd00016==0 , cpd00024==0 , cpd00040==0 , cpd00076==0 , cpd00077==0 , cpd00098==0 , cpd00104==0 , cpd00118==0 , cpd00121==0 , cpd00122==0 , cpd00136==0 , cpd00158==0 , cpd00160==0 , cpd00162==0 , cpd00164==0 , cpd00166==0 , cpd00179==0 , cpd00208==0 , cpd00210==0 , cpd00215==0 , cpd00218==0 , cpd00220==0 , cpd00222==0 , cpd00226==0 , cpd00228==0 , cpd00240==0 , cpd00245==0 , cpd00246==0 , cpd00247==0 , cpd00248==0 , cpd00250==0 , cpd00263==0 , cpd00280==0 , cpd00281==0 , cpd00300==0 , cpd00305==0 , cpd00314==0 , cpd00324==0 , cpd00328==0 , cpd00332==0 , cpd00338==0 , cpd00339==0 , cpd00359==0 , cpd00361==0 , cpd00382==0 , cpd00393==0 , cpd00396==0 , cpd00419==0 , cpd00420==0 , cpd00441==0 , cpd00443==0 , cpd00453==0 , cpd00456==0 , cpd00479==0 , cpd00504==0 , cpd00528==0 , cpd00536==0 , cpd00540==0 , cpd00541==0 , cpd00585==0 , cpd00588==0 , cpd00599==0 , cpd00601==0 , cpd00636==0 , cpd00644==0 , cpd00646==0 , cpd00666==0 , cpd00680==0 , cpd00794==0 , cpd00813==0 , cpd00992==0 , cpd00996==0 , cpd01020==0 , cpd01034==0 , cpd01059==0 , cpd01089==0 , cpd01092==0 , cpd01208==0 , cpd01223==0 , cpd01313==0 , cpd01401==0 , cpd01415==0 , cpd01511==0 , cpd01530==0 , cpd01574==0 , cpd01582==0 , cpd01583==0 , cpd01679==0 , cpd01685==0 , cpd01703==0 , cpd01706==0 , cpd01822==0 , cpd01826==0 , cpd01835==0 , cpd01947==0 , cpd02128==0 , cpd02197==0 , cpd02246==0 , cpd02301==0 , cpd03026==0 , cpd03185==0 , cpd03292==0 , cpd03533==0 , cpd03644==0 , cpd03846==0 , cpd03959==0 , cpd03962==0 , cpd03977==0 , cpd03981==0 , cpd03998==0 , cpd04018==0 , cpd04019==0 , cpd04073==0 , cpd04078==0 , cpd04135==0 , cpd04145==0 , cpd04148==0 , cpd04184==0 , cpd04307==0 , cpd04309==0 , cpd04358==0 , cpd04367==0 , cpd04371==0 , cpd04373==0 , cpd04375==0 , cpd04446==0 , cpd04450==0 , cpd04623==0 , cpd04733==0 , cpd04864==0 , cpd05178==0 , cpd06523==0 , cpd07326==0 , cpd07339==0 , cpd07719==0 , cpd07891==0 , cpd07905==0 , cpd08020==0 , cpd08027==0 , cpd08053==0 , cpd08191==0 , cpd08276==0 , cpd09057==0 , cpd09225==0 , cpd09244==0 , cpd09248==0 , cpd09320==0 , cpd09918==0 , cpd10023==0 , cpd10025==0 , cpd10027==0 , cpd10095==0 , cpd10107==0 , cpd10149==0 , cpd10155==0 , cpd10232==0 , cpd10393==0 , cpd10908==0 , cpd11145==0 , cpd11181==0 , cpd11575==0 , cpd11594==0 , cpd11601==0 , cpd11602==0 , cpd11624==0 , cpd11645==0 , cpd11656==0 , cpd11657==0 , cpd11658==0 , cpd11664==0 , cpd11683==0 , cpd11688==0 , cpd11732==0 , cpd11746==0 , cpd11862==0 , cpd11879==0 , cpd11904==0 , cpd11949==0 , cpd12183==0 , cpd12773==0 , cpd12836==0 , cpd12859==0 , cpd13334==0 , cpd13391==0 , cpd13392==0 , cpd14465==0 , cpd14795==0 , cpd14913==0 , cpd14921==0 , cpd15142==0 , cpd15608 ==0)
#Only then can we remove those compound columns from consideration.
print('...then removing those compounds from consideration...')
filtered_media <- select(filtered_media, -c( cpd00010, cpd00015, cpd00016, cpd00024, cpd00040, cpd00076, cpd00077, cpd00098, cpd00104, cpd00118, cpd00121, cpd00122, cpd00136, cpd00158, cpd00160, cpd00162, cpd00164, cpd00166, cpd00179, cpd00208, cpd00210, cpd00215, cpd00218, cpd00220, cpd00222, cpd00226, cpd00228, cpd00240, cpd00245, cpd00246, cpd00247, cpd00248, cpd00250, cpd00263, cpd00280, cpd00281, cpd00300, cpd00305, cpd00314, cpd00324, cpd00328, cpd00332, cpd00338, cpd00339, cpd00359, cpd00361, cpd00382, cpd00393, cpd00396, cpd00419, cpd00420, cpd00441, cpd00443, cpd00453, cpd00456, cpd00479, cpd00504, cpd00528, cpd00536, cpd00540, cpd00541, cpd00585, cpd00588, cpd00599, cpd00601, cpd00636, cpd00644, cpd00646, cpd00666, cpd00680, cpd00794, cpd00813, cpd00992, cpd00996, cpd01020, cpd01034, cpd01059, cpd01089, cpd01092, cpd01208, cpd01223, cpd01313, cpd01401, cpd01415, cpd01511, cpd01530, cpd01574, cpd01582, cpd01583, cpd01679, cpd01685, cpd01703, cpd01706, cpd01822, cpd01826, cpd01835, cpd01947, cpd02128, cpd02197, cpd02246, cpd02301, cpd03026, cpd03185, cpd03292, cpd03533, cpd03644, cpd03846, cpd03959, cpd03962, cpd03977, cpd03981, cpd03998, cpd04018, cpd04019, cpd04073, cpd04078, cpd04135, cpd04145, cpd04148, cpd04184, cpd04307, cpd04309, cpd04358, cpd04367, cpd04371, cpd04373, cpd04375, cpd04446, cpd04450, cpd04623, cpd04733, cpd04864, cpd05178, cpd06523, cpd07326, cpd07339, cpd07719, cpd07891, cpd07905, cpd08020, cpd08027, cpd08053, cpd08191, cpd08276, cpd09057, cpd09225, cpd09244, cpd09248, cpd09320, cpd09918, cpd10023, cpd10025, cpd10027, cpd10095, cpd10107, cpd10149, cpd10155, cpd10232, cpd10393, cpd10908, cpd11145, cpd11181, cpd11575, cpd11594, cpd11601, cpd11602, cpd11624, cpd11645, cpd11656, cpd11657, cpd11658, cpd11664, cpd11683, cpd11688, cpd11732, cpd11746, cpd11862, cpd11879, cpd11904, cpd11949, cpd12183, cpd12773, cpd12836, cpd12859, cpd13334, cpd13391, cpd13392, cpd14465, cpd14795, cpd14913, cpd14921, cpd15142, cpd15608 ))
print('done!')

#...then label compounds found in CHNOSZ, but which are not present in the media.
print('Removing compounds not actually used in any media...')
cpds <- data.frame("nonzero" = colSums((filtered_media)[-c(1,2)])==0)
cpds['cpd'] <- rownames(cpds)
cpds <- filter(cpds, nonzero==TRUE)$cpd
cpds_not_used_2 <- paste('filtered_media <- select(filtered_media, -c(', paste(cpds, collapse=', '), '))')
filtered_media <- select(filtered_media, -c( cpd00007, cpd00011, cpd00012, cpd00021, cpd00033, cpd00035, cpd00039, cpd00047, cpd00060, cpd00066, cpd00067, cpd00075, cpd00106, cpd00107, cpd00128, cpd00129, cpd00133, cpd00139, cpd00141, cpd00154, cpd00159, cpd00161, cpd00178, cpd00184, cpd00204, cpd00207, cpd00211, cpd00224, cpd00308, cpd00322, cpd00367, cpd00379, cpd00386, cpd00411, cpd00450, cpd00547, cpd00552, cpd00572, cpd00618, cpd00637, cpd00797, cpd00998, cpd01007, cpd01024, cpd01041, cpd01042, cpd01048, cpd01079, cpd01080, cpd01087, cpd01112, cpd01269, cpd01391, cpd01642, cpd01741, cpd01834, cpd03387, cpd03396, cpd03836, cpd03994, cpd04166, cpd04877, cpd05267, cpd09249, cpd09400, cpd09693, cpd09695, cpd11161, cpd11640, cpd15574 ))
print('done!')


### Replace compound numbers with CHNOSZ names
##### Keep in mind that SEED names sometimes have stereochemistry embedded, but CHNOSZ doesn't distinguish enantiomers.

filtered_media <- filter(media_usable,  cpd00010==0 , cpd00015==0 , cpd00016==0 , cpd00024==0 , cpd00040==0 , cpd00076==0 , cpd00077==0 , cpd00098==0 , cpd00104==0 , cpd00118==0 , cpd00121==0 , cpd00122==0 , cpd00136==0 , cpd00158==0 , cpd00160==0 , cpd00162==0 , cpd00164==0 , cpd00166==0 , cpd00179==0 , cpd00208==0 , cpd00210==0 , cpd00215==0 , cpd00218==0 , cpd00220==0 , cpd00222==0 , cpd00226==0 , cpd00228==0 , cpd00240==0 , cpd00245==0 , cpd00246==0 , cpd00247==0 , cpd00248==0 , cpd00250==0 , cpd00263==0 , cpd00280==0 , cpd00281==0 , cpd00300==0 , cpd00305==0 , cpd00314==0 , cpd00324==0 , cpd00328==0 , cpd00332==0 , cpd00338==0 , cpd00339==0 , cpd00359==0 , cpd00361==0 , cpd00382==0 , cpd00393==0 , cpd00396==0 , cpd00419==0 , cpd00420==0 , cpd00441==0 , cpd00443==0 , cpd00453==0 , cpd00456==0 , cpd00479==0 , cpd00504==0 , cpd00528==0 , cpd00536==0 , cpd00540==0 , cpd00541==0 , cpd00585==0 , cpd00588==0 , cpd00599==0 , cpd00601==0 , cpd00636==0 , cpd00644==0 , cpd00646==0 , cpd00666==0 , cpd00680==0 , cpd00794==0 , cpd00813==0 , cpd00992==0 , cpd00996==0 , cpd01020==0 , cpd01034==0 , cpd01059==0 , cpd01089==0 , cpd01092==0 , cpd01208==0 , cpd01223==0 , cpd01313==0 , cpd01401==0 , cpd01415==0 , cpd01511==0 , cpd01530==0 , cpd01574==0 , cpd01582==0 , cpd01583==0 , cpd01679==0 , cpd01685==0 , cpd01703==0 , cpd01706==0 , cpd01822==0 , cpd01826==0 , cpd01835==0 , cpd01947==0 , cpd02128==0 , cpd02197==0 , cpd02246==0 , cpd02301==0 , cpd03026==0 , cpd03185==0 , cpd03292==0 , cpd03533==0 , cpd03644==0 , cpd03846==0 , cpd03959==0 , cpd03962==0 , cpd03977==0 , cpd03981==0 , cpd03998==0 , cpd04018==0 , cpd04019==0 , cpd04073==0 , cpd04078==0 , cpd04135==0 , cpd04145==0 , cpd04148==0 , cpd04184==0 , cpd04307==0 , cpd04309==0 , cpd04358==0 , cpd04367==0 , cpd04371==0 , cpd04373==0 , cpd04375==0 , cpd04446==0 , cpd04450==0 , cpd04623==0 , cpd04733==0 , cpd04864==0 , cpd05178==0 , cpd06523==0 , cpd07326==0 , cpd07339==0 , cpd07719==0 , cpd07891==0 , cpd07905==0 , cpd08020==0 , cpd08027==0 , cpd08053==0 , cpd08191==0 , cpd08276==0 , cpd09057==0 , cpd09225==0 , cpd09244==0 , cpd09248==0 , cpd09320==0 , cpd09918==0 , cpd10023==0 , cpd10025==0 , cpd10027==0 , cpd10095==0 , cpd10107==0 , cpd10149==0 , cpd10155==0 , cpd10232==0 , cpd10393==0 , cpd10908==0 , cpd11145==0 , cpd11181==0 , cpd11575==0 , cpd11594==0 , cpd11601==0 , cpd11602==0 , cpd11624==0 , cpd11645==0 , cpd11656==0 , cpd11657==0 , cpd11658==0 , cpd11664==0 , cpd11683==0 , cpd11688==0 , cpd11732==0 , cpd11746==0 , cpd11862==0 , cpd11879==0 , cpd11904==0 , cpd11949==0 , cpd12183==0 , cpd12773==0 , cpd12836==0 , cpd12859==0 , cpd13334==0 , cpd13391==0 , cpd13392==0 , cpd14465==0 , cpd14795==0 , cpd14913==0 , cpd14921==0 , cpd15142==0 , cpd15608 ==0)
filtered_media <- select(filtered_media, -c( cpd00010, cpd00015, cpd00016, cpd00024, cpd00040, cpd00076, cpd00077, cpd00098, cpd00104, cpd00118, cpd00121, cpd00122, cpd00136, cpd00158, cpd00160, cpd00162, cpd00164, cpd00166, cpd00179, cpd00208, cpd00210, cpd00215, cpd00218, cpd00220, cpd00222, cpd00226, cpd00228, cpd00240, cpd00245, cpd00246, cpd00247, cpd00248, cpd00250, cpd00263, cpd00280, cpd00281, cpd00300, cpd00305, cpd00314, cpd00324, cpd00328, cpd00332, cpd00338, cpd00339, cpd00359, cpd00361, cpd00382, cpd00393, cpd00396, cpd00419, cpd00420, cpd00441, cpd00443, cpd00453, cpd00456, cpd00479, cpd00504, cpd00528, cpd00536, cpd00540, cpd00541, cpd00585, cpd00588, cpd00599, cpd00601, cpd00636, cpd00644, cpd00646, cpd00666, cpd00680, cpd00794, cpd00813, cpd00992, cpd00996, cpd01020, cpd01034, cpd01059, cpd01089, cpd01092, cpd01208, cpd01223, cpd01313, cpd01401, cpd01415, cpd01511, cpd01530, cpd01574, cpd01582, cpd01583, cpd01679, cpd01685, cpd01703, cpd01706, cpd01822, cpd01826, cpd01835, cpd01947, cpd02128, cpd02197, cpd02246, cpd02301, cpd03026, cpd03185, cpd03292, cpd03533, cpd03644, cpd03846, cpd03959, cpd03962, cpd03977, cpd03981, cpd03998, cpd04018, cpd04019, cpd04073, cpd04078, cpd04135, cpd04145, cpd04148, cpd04184, cpd04307, cpd04309, cpd04358, cpd04367, cpd04371, cpd04373, cpd04375, cpd04446, cpd04450, cpd04623, cpd04733, cpd04864, cpd05178, cpd06523, cpd07326, cpd07339, cpd07719, cpd07891, cpd07905, cpd08020, cpd08027, cpd08053, cpd08191, cpd08276, cpd09057, cpd09225, cpd09244, cpd09248, cpd09320, cpd09918, cpd10023, cpd10025, cpd10027, cpd10095, cpd10107, cpd10149, cpd10155, cpd10232, cpd10393, cpd10908, cpd11145, cpd11181, cpd11575, cpd11594, cpd11601, cpd11602, cpd11624, cpd11645, cpd11656, cpd11657, cpd11658, cpd11664, cpd11683, cpd11688, cpd11732, cpd11746, cpd11862, cpd11879, cpd11904, cpd11949, cpd12183, cpd12773, cpd12836, cpd12859, cpd13334, cpd13391, cpd13392, cpd14465, cpd14795, cpd14913, cpd14921, cpd15142, cpd15608 ))
filtered_media <- select(filtered_media, -c( cpd00007, cpd00011, cpd00012, cpd00021, cpd00033, cpd00035, cpd00039, cpd00047, cpd00060, cpd00066, cpd00067, cpd00075, cpd00106, cpd00107, cpd00128, cpd00129, cpd00133, cpd00139, cpd00141, cpd00154, cpd00159, cpd00161, cpd00178, cpd00184, cpd00204, cpd00207, cpd00211, cpd00224, cpd00308, cpd00322, cpd00367, cpd00379, cpd00386, cpd00411, cpd00450, cpd00547, cpd00552, cpd00572, cpd00618, cpd00637, cpd00797, cpd00998, cpd01007, cpd01024, cpd01041, cpd01042, cpd01048, cpd01079, cpd01080, cpd01087, cpd01112, cpd01269, cpd01391, cpd01642, cpd01741, cpd01834, cpd03387, cpd03396, cpd03836, cpd03994, cpd04166, cpd04877, cpd05267, cpd09249, cpd09400, cpd09693, cpd09695, cpd11161, cpd11640, cpd15574 ))

`%notin%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))
i=3

#Double check the stereochemistry - combine both forms into one column
while (i<=ncol(filtered_media)){
    colname=names(filtered_media)[i]
    chnoszID <- paste(filter(CHNOSZ_conv, id==colname)$CHNOSZ_id)
    SEEDid <- filter(CHNOSZ_conv, id==colname)
    if (chnoszID %notin% colnames(filtered_media)) {
        names(filtered_media)[i] <- chnoszID
    }
    else {  #if applicable, combines L and D enantiomers, or other species possibly "identical" in CHNOSZ
        filtered_media[chnoszID] <- filtered_media[chnoszID] + filtered_media[i]
    }
    i=i+1
}
#Remove relic columns from (excess) stereochemical forms - they've already been combined in CHNOSZ-named fields
reduced_media <- filtered_media %>% select(-starts_with("cpd"))
write.csv(reduced_media, "reduced_media.csv")
print(filter(data.frame(summary(reduced_media)), grepl("Max",Freq))) #Show max molality - i.e. how important is it in at least one case?

speciesname <- (names(reduced_media))[-c(1,2)]
species <- data.frame(speciesname)
queries <- data.frame(check.names=TRUE, name=character(), abbrv=character(), formula=character(), state=character())
for (i in species$speciesname){
    query1 <- filter(CHNOSZ$thermo$obigt, str_detect(name, paste(i))==TRUE)[c(1,2,3,4)] #$#, str_detect(abbrv, 'NAD')==TRUE)'
    query2 <- filter(CHNOSZ$thermo$obigt, str_detect(name, paste(i))==FALSE, str_detect(abbrv, paste(i))==TRUE)[c(1,2,3,4)]
	poss_relat <- str_to_title(paste(substr(i[1], start=1, stop=3)))
	#cap1 <- str_to_title(paste(i[1:3]))
    query3 <- filter(CHNOSZ$thermo$obigt, str_detect(name, paste(i))==FALSE, str_detect(abbrv, paste(i))==FALSE, str_detect(abbrv, paste(poss_relat))==TRUE)[c(1,2,3,4)]
    queries <- rbind(queries, query1, query2, query3)
    queries <- queries[!duplicated(queries), ]
}

write.csv(queries, "candidate_compounds.csv")
