# Initial script created by Arthur Elmes, with significant contributions from
# A.J. Shatz, based on methods of Jimenez-Munoz et al.

# This thing needs to be really reworked, ideally abandoning arcpy in favor of
# gdal or rasterio, and better yet made into a module so the functions can be
# invoked as-needed, rather than as a big waterfall

# import stuff
import os
import shutil
import arcpy
from arcpy import env
from arcpy.sa import *
# import win32com.client
import math
import glob
from glob import glob
import csv
import numpy as np
import datetime

# print processing start time
print "Processing started: ",
print datetime.datetime.now()
startTime = datetime.datetime.now()

# Check out spatial analyst license (why?)
arcpy.CheckOutExtension("Spatial")

# set directories
# the 'test' directories are all identical, but are on D instead of G

working_folder = "C:\\Data\\Landsat"  # type: str
refl_dl_folder = os.path.join(working_folder, "OriginalDLFolder\\Reflectance")
lvl1_dl_folder = os.path.join(working_folder, "OriginalDLFolder\\LVL1")
output_folder_root = os.path.join(working_folder, "Output")
# DELETE?? refl_output_folder = os.path.join(working_folder, "Reflectance")

# arcpy environment settings (blerg)
env.overwriteOutput = True
env.workspace = output_folder_root
# DELETE??? tile_folder = env.workspace

# list of path/row combos and sensors
path_row_list = ["012031", "013030", "013031"]  # might add to these if necessary
sensor_list = ["LT05", "LE07", "LC08"]

# this is the list of rasters that, for whatever reason, fail to be correctly processed
error_list_location = "C:\\Data\\Landsat\\error_list.txt"

# clear any rasters from previous run
error_list = open(error_list_location, "w")
error_list.close()


# Loop through 'original' folder containing all LC8, LE7, and LT5 scenes from the CDR archive
# Sort the raw download data into the appropriate folders, organized first by tile
# (path/row) then by sensor (LT5, LE7, or LC8)
def sort_images(lvl1_download_folder, refl_download_folder, refl_out_folder, pr_list, sns_list):
    # First, create any pathrow folders needed in the output directory, based on the pr and sensor lists provided above
    for pr in pr_list:
        for sensor in sns_list:
            out_dir = os.path.join(refl_out_folder, pr, sensor)
            if not os.path.isdir(out_dir):
                os.mkdir(out_dir)

    # Then, move the files into their correct folder
    for img_folder in os.listdir(refl_download_folder):
        tile_name = str(img_folder)
        for pr in pr_list:
            if pr in tile_name:
                for sns in sns_list:
                    if sns in tile_name:
                        print tile_name, "is", sns, "and is in", pr
                        src_fldr = os.path.join(refl_download_folder, tile_name)
                        dest_fldr = os.path.join(refl_out_folder, pr, sns, tile_name)
                        if not os.path.isdir(dest_fldr):
                            shutil.copytree(src_fldr, dest_fldr)
                        find_thermal(lvl1_download_folder, dest_fldr)


# Find the raw thermal band for each image
# Loop through all CDR images, and for each one, loop through the LVL1 images, stopping at the matching image,
# based on name
# Copy that image and paste it in the CDR folder


def find_thermal(lvl1_download_folder, refl_dest_folder):
    # parse the file name
    # THIS MAY NEED TO CHANGE, since apparently these conventions change occasionally
    found = False
    sensor = refl_dest_folder.split("\\")[6][:4]
    pathrow = refl_dest_folder.split("\\")[6][4:10]
    for root, dirs, files in os.walk(lvl1_download_folder):
        for dir in dirs:
            if pathrow in dir and sensor in dir:
                found = True
                imgs = os.listdir(os.path.join(root, dir))
                for img in imgs:
                    if "LC08" in img:
                        if "B10" in img:
                            shutil.copy2(os.path.join(root, dir, img), refl_dest_folder)
                        elif "MTL" in dir:
                            shutil.copy2(os.path.join(root, dir, img), refl_dest_folder)
                    elif "LT05" in sensor or "LE07" in sensor:
                        if "B6" in img or "VCID" in img:
                            shutil.copy2(os.path.join(root, dir, img), refl_dest_folder)
                        elif "MTL" in img:
                            shutil.copy2(os.path.join(root, dir, img), refl_dest_folder)
    if not found:
        error_list = open(error_list_location, "a")
        error_list.write("No thermal image found for: " + str(refl_dest_folder))
        error_list.close()


sort_images(lvl1_dl_folder, refl_dl_folder, output_folder_root, path_row_list, sensor_list)


# iterate through each band in the folder, and assign raster objects to the cloud mask,
# NIR, and red bands
def set_masks():
    try:
        for tile in os.listdir(tileFolderRoot):
            for sensor in os.listdir(os.path.join(tileFolderRoot, tile)):
                for image in os.listdir(os.path.join(tileFolderRoot, tile, sensor)):
                    env.workspace = os.path.join(tileFolderRoot, tile, sensor, image)
                    imageFolder = os.path.join(tileFolderRoot, tile, sensor, image)
                    bandList = arcpy.ListRasters("*", "TIF")
                    for band in bandList:
                        if "LC8" in band:
                            if "band4" in band:
                                redBand = arcpy.Raster(imageFolder + "\\" + band)
                            elif "band5" in band:
                                nirBand = arcpy.Raster(imageFolder + "\\" + band)
                            elif "cfmask.tif" in band:
                                cloudMaskOrig = arcpy.Raster(imageFolder + "\\" + band)
                        elif "LE7" in band:
                            if "band3" in band:
                                redBand = arcpy.Raster(imageFolder + "\\" + band)
                            elif "band4" in band:
                                nirBand = arcpy.Raster(imageFolder + "\\" + band)
                            elif "cfmask.tif" in band:
                                cloudMaskOrig = arcpy.Raster(imageFolder + "\\" + band)
                        elif "LT5" in band:
                            if "band3" in band:
                                redBand = arcpy.Raster(imageFolder + "\\" + band)
                            elif "band4" in band:
                                nirBand = arcpy.Raster(imageFolder + "\\" + band)
                            elif "cfmask.tif" in band:
                                cloudMaskOrig = arcpy.Raster(imageFolder + "\\" + band)
                    redBandFloat = arcpy.sa.Float(redBand)
                    nirBandFloat = arcpy.sa.Float(nirBand)
                    print "Calculating NDVI for image: " + image + " : " + redBand.name + " and " + nirBand.name
                    try:
                        ndviImage = Float(nirBand - redBand)/(nirBand + redBand)
                        print "image folder and name: ",
                        print imageFolder + "\\" + bandList[0][:16] + "_NDVI.tif"
                        ndviImage.save(imageFolder + "\\" + bandList[0][:16] + "_NDVI.tif")
                    except:
                        error_list = open(error_listLocation, "a")
                        error_list.write("NDVI calculation fail: " + str(nirBand))
                        error_list.close()

                    # take the mask and create a binary mask to remove all snow, cloud, cloud shadow, and water.
                    # Store the remap values in a variable called remapKey
                    # remapKey = RemapValue([[0,1],[1,"NODATA"],[2,"NODATA"],[3,"NODATA"],[4,"NODATA"]])
                    try:
                        remapKey = RemapValue([[2,1],[3,1],[4,1]])
                        cloudMaskBin = Reclassify(cloudMaskOrig, "Value", remapKey)
                        print "Creating cloudmask: ",
                        print cloudMaskBin
                        try:
                            cloudMaskBin.save(cloudMaskOrig.path + cloudMaskOrig.name[:-4] + "_binary.tif")
                        except:
                            print "Can't create cloudmask for: ",
                            print str(cloudMaskBin)
                    except:
                        error_list = open(error_listLocation, "a")
                        error_list.write("CloudMask fail: " + str(cloudMaskBin))
                        error_list.close()
    except:
        error_list = open(error_listLocation, "a")
        error_list.write("CloudMask fail: " + str(tile))
        error_list.close()
        print "Could not process: ",
        print str(tile)


# this deletes any previously made rasters, just so no duplicates are made
def delete_rasters(currentFolder, currentFileList):
    env.workspace = currentFolder
    tempRasterList = arcpy.ListRasters("*","TIF")
    for rasters in tempRasterList:
        if "lSen" in rasters or "tSen" in rasters or "gamma" in rasters or "delta" in rasters or "emiss" in rasters or "lst" in rasters or "NDVI" in rasters:
            if "Masked" not in rasters:
                arcpy.Delete_management(rasters)
    del tempRasterList


# these are the modeled radiosounding values from the JMS series of papers
radiosoundDict = {'LT5':[[0.08735,-0.09553,1.10188],[-0.69188,-0.58185,-0.29887],[-0.03724,1.53065,-0.45476]],\
                  'LE7':[[0.07593,-0.07132,1.08565],[-0.61438,-0.70916,-0.19379],[-0.02892,1.46051,-0.43199]],\
                  'LC8':[[0.04019,0.02916,1.01523],[-0.38333,-1.50294,0.20324],[0.00918,1.36072,-0.27514]]}

# This function calculates at sensor radiance by multiplying the image by the gain and adding the bais.
# These items are retrieved from the mtl metadata file.
def lSen(thermalImageFolder, thermalImage, metaFile):
    imageName = str(thermalImage)
    metaFilePath = ""
    gain = ""
    bias = ""
    lSenRaster = ""
    if len(imageName) > 0:
        if "LT5" in imageName:
            if len(metaFile) > 0:
                metadata = open(metaFile, 'r')
                for line in metadata:
                    if "RADIANCE_MULT_BAND_6" in line:
                        gain = float(str.split(line)[2])
                    elif "RADIANCE_ADD_BAND_6" in line:
                        bias = float(str.split(line)[2])
        elif "LE7" in imageName:
            if len(metaFile) > 0:
                metadata = open(metaFile, 'r')
                for line in metadata:
                    if "RADIANCE_MULT_BAND_6_VCID_1" in line:
                        gain = float(str.split(line)[2])
                    elif "RADIANCE_ADD_BAND_6_VCID_1" in line:
                        bias = float(str.split(line)[2])
        elif "LC8" in imageName:
            if len(metaFile) > 0:
                metadata = open(metaFile, 'r')
                for line in metadata:
                    if "RADIANCE_MULT_BAND_10" in line:
                        gain = float(str.split(line)[2])
                    elif "RADIANCE_ADD_BAND_10" in line:
                        bias = float(str.split(line)[2])
        try:
            lSenRaster = bias + Raster(thermalImageFolder + "\\" + imageName) * gain
            # print thermalImageFolder + "\\" + imageName[:-4] + "_lSen.TIF"
            lSenRaster.save(thermalImageFolder + "\\" + imageName[:-4] + "_lSen.TIF")
            return lSenRaster.path + lSenRaster.name
        except:
            error_list = open(error_listLocation, "a")
            error_list.write("lSen: " + str(thermalImage))
            error_list.close()

    return lSenRaster


def tSen(lSenRasterFolder, lSenRaster):
    try:
        env.workspace = lSenRasterFolder
        radRasterName = str(lSenRaster)
        k1 = ""
        k2 = ""
        tSenRaster = ""
        if len(radRasterName) > 0 and radRasterName <> "None":
            if "LT5" in radRasterName:
                k1 = 607.76
                k2 = 1260.56
            elif "LE7" in radRasterName:
                k1 = 666.09
                k2 = 1282.71
            elif "LC8" in radRasterName:
                k1 = 774.89
                k2 = 1321.08
            tSenRaster = k2 / arcpy.sa.Ln(k1 / arcpy.Raster(radRasterName) + 1)
            tSenRaster.save(str(lSenRaster)[:-8] + "_tSen.TIF")
            #print str(lSenRaster)[:-8] + "_tSen.TIF"
            return tSenRaster.path + tSenRaster.name
            del tSenRaster
            del radRastername
    except:
        print arcpy.AddError(arcpy.GetMessages(2))
        print "Could not calculate tSen for: ",
        print lSenRaster
        error_list = open(error_listLocation, "a")
        error_list.write("tSen calculation fail: " + str(lSenRaster))
        error_list.close()

def gammaDelta(currentFolder, tSenRaster, lSenRaster):
    try:
        env.workspace = currentFolder
        tSenRasterName = str(tSenRaster)
        lamb = 0.0
        if len(str(tSenRasterName)) > 0 and "None" not in tSenRasterName:
            # these lambda coefficients are taken from Jimenez et al 2003, 2009, and 2014
            # check the lambdas for units -- AJ's were two orders smaller
            # I actually made these two orders of magnitude bigger than JMS
            if "LT5" in tSenRasterName:
                lamb = 125600
            elif "LE7" in tSenRasterName:
                lamb = 127700
            elif "LC8" in tSenRasterName:
                lamb = 132400
            # calculate gamma and delta, as defined by Jimenez et al 2003, 2009, and 2014
            gammaRaster = arcpy.sa.Power(tSenRaster, 2) / lamb * Raster(lSenRaster)
            deltaRaster = Raster(tSenRaster) - arcpy.sa.Power(tSenRaster, 2) / lamb
            gammaRaster.save(str(tSenRasterName)[:-9] + "gamma.TIF")
            deltaRaster.save(str(tSenRasterName)[:-9] + "delta.TIF")
            return gammaRaster.name, deltaRaster.name
            del gammaRaster
            del deltaRaster
    except:
        print arcpy.AddError(arcpy.GetMessages(2))
        print "Could not calculate Gamma and Delta for: ",
        print lSenRaster
        error_list = open(error_listLocation, "a")
        error_list.write("gamma or delta calculation fail: " + str(lSenRaster))
        error_list.close()

def emissivity(currentFolder, ndvi):
    # print "Calculating emissivity for: ",
    # print ndvi
    try:
        env.workspace = currentFolder
        ndviRasterName = str(ndvi)
        emisRaster = ""
        emisRasterReclass = ""
        propVegetation = ""
        maskBand = ""
        bandList = arcpy.ListRasters("*", "TIF")
        for band in bandList:
            if "_binary.tif" in band:
                maskBand = arcpy.Raster(band)
        if len(ndviRasterName) > 0 and "None" not in ndviRasterName:
            # reclasfiy the NDVI image to emissivity based on Sobrino et al 2004 thresholds of
            # <0.2 for bare soil, >0.5 for fully vegetated, and 0.2<NDVI<0.5 being a mixed pixel
            # LSE values also came from Sobrino et al 2004
            # This needs to be reworked -- the MXD on the desktop has the correct steps
            # Basically the emissivity raster will be created by raster adding three rasters
            # First CON the soil and veg emissivity rasters: NDVI < 0.2 = emis 0.973; NDVI> 0.5 = emis 0.99
            # for all other cases, reclass to 0
            # Then CON the mixed soil/veg emis by including a conditoinal that only calculates new emis values
            # when NDVI are within the range, otherwise reclass to 0
            # This should yeild three interlocking rasters, which can just be added together for final emissivity.
            vegEmisRst = arcpy.sa.Con(ndviRasterName, 0.99, 0, "VALUE >= 0.5")
            soilEmisRst = arcpy.sa.Con(ndviRasterName, 0.973, 0, "VALUE <= 0.2")
            mixedEmisRst = arcpy.sa.Con(ndviRasterName, 0.004 * arcpy.sa.Power((Raster(ndviRasterName) - 0.2/0.5-0.2),2) + 0.986, 0, "VALUE > 0.2 AND VALUE < 0.5")
            combinedEmisRst = mixedEmisRst + soilEmisRst + vegEmisRst
            # combinedEmisRst = combinedEmisRst * maskBand
            combinedEmisRst = arcpy.sa.Con(maskBand, combinedEmisRst, -999, "value = 0")
            print str(ndviRasterName[:-14] + "_emissivity_SA.TIF")
            combinedEmisRst.save(str(ndviRasterName[:-14] + "_emissivity_SA.TIF"))
        return combinedEmisRst.path + combinedEmisRst.name

    except:
        print arcpy.AddError(arcpy.GetMessages(2))
        print "Could not calculate emissivity for: ",
        print currentFolder
        error_list = open(error_listLocation, "a")
        error_list.write("emissivity calculation fail, substituted placeHolderRst: " + str(ndvi))
        error_list.close()
        placeHolderRst = arcpy.Raster("G:\\Data\\Dissertation\\LandsatDataWorcester\\PlaceHolderRaster\\PlaceHolder_AllNoData.tif")
        # placeHolderRst.save(os.path.join(currentFolder, currentFolder[-33:-17] + "_emissivity_SA.TIF"))
        imageNameList = currentFolder.split(os.sep)
        imageName = imageNameList[len(imageNameList) - 1]
        placeHolderRstSavePath = os.path.join(currentFolder, imageName + "_emissivity_SA.TIF")
        placeHolderRst.save(placeHolderRstSavePath)
        print "Placeholder saved: ",
        print placeHolderRstSavePath
        return placeHolderRstSavePath

def w(currentImage):
    # retrieve the relative humidity value fom the worksheet defined below,
    # which contains the metadata for calculation and data sources
    imageDate = str(currentImage)[9:16]
    relHum = 0.0
    with open(r"D:\sync\gis\Dissertation\scripts\Thermal\WeatherDataWorcester.csv", "rb") as csvfile:
       dataReader = csv.reader(csvfile)
       for row in dataReader:
           if imageDate in row:
               relHum = float(row[3]) * 0.01
    return relHum
    del relHum

def psi123(radiosoundDict, w, thermalBand):
    psi = []
    if len(str(thermalBand)) > 0:
        # coefs = [[1],[w],[w**2]]
        coefs = [[w**2],[w],[1]]
        coefMatrix = np.matrix(coefs)
        sens = str(thermalBand)[0:3]
        radiosoundValueMatrix = np.matrix(radiosoundDict[sens])
        psi = radiosoundValueMatrix * coefMatrix
        return psi
        del psi

def lst(currentFolder, lSenRaster, gammaDeltaList, emisRaster, psiValues):
    try:
        env.workspace = currentFolder
        lstRaster = arcpy.sa.Con(emisRaster, arcpy.sa.Float(Raster(gammaDeltaList[0]) * (1/Raster(emisRaster)\
                    * (float(psiValues[0]) * (Raster(lSenRaster) * 0.1) + float(psiValues[1])) + float(psiValues[2])) + Raster(gammaDeltaList[1])), -999, "value > -999")
        # lstRaster.save("D:\\Data\\Dissertation\\LandsatDataWorcester\\LST_Run3\\" + str(lSenRaster[-36:-8]) + "lst.TIF")
        # lstRaster.save("D:\\Data\\Dissertation\\LandsatDataWorcester\\LST_Run3\\" + lSenRaster.name[:-8] + "lst.TIF")

        # get the correct image name from the current folder name, minus the CDR-specific gobldegook.
        imageNameList = currentFolder.split(os.sep)
        imageName1 = imageNameList[len(imageNameList) - 1]
        imageName2 = imageName1.split("-")
        if len(imageName2) > 0:
            imageName3 = imageName2[0]
        else:
            print "Woops"
        print imageName3
        lstRasterSavePath = os.path.join("D:\\Data\\Dissertation\\LandsatDataWorcester\\LST_Run3\\RawLSTSeries", imageName3 + "_lst.TIF")
        print "Saving LST to: ",
        print lstRasterSavePath
        lstRaster.save(lstRasterSavePath)
        # lstRaster.save(str(lSenRaster + "_lst.TIF"))
        del lstRaster, gammaDeltaList, emisRaster, psiValues, lstRasterSavePath, imageNameList, imageName1, imageName2, imageName3

    except:
        print "Could not calculate lst for: ",
        print lSenRaster
        error_list = open(error_listLocation, "a")
        error_list.write("LST calculation fail within module: " + str(lSenRaster))
        error_list.close()


# Mask the lst and NDVI images
def maskBands(currentFolder, currentFileList):
    env.workspace = currentFolder
    bandList = arcpy.ListRasters("*", "TIF")
    maskBand = ""
    for band in bandList:
        if "_binary.tif" in band:
            maskBand = arcpy.Raster(band)
    for bandToMask in bandList:
        if "lst.TIF" in bandToMask:
            if "Masked" not in bandToMask:
                print "Masking image: "
                toMask = arcpy.Raster(bandToMask)
                print  bandToMask
                maskedBand = arcpy.sa.Con(maskBand, toMask, "-999", "Value = 0")
                # maskedBand = toMask * maskBand
                maskedBand.save(toMask.name[:-4]  + "_Masked.tif")
                del toMask
        elif "NDVI" in bandToMask and "emissivity" not in bandToMask:
            if "Masked" not in bandToMask:
                print "Masking image: "
                toMask = arcpy.Raster(bandToMask)
                print  bandToMask
                # maskedBand = toMask * maskBand
                maskedBand = arcpy.sa.Con(maskBand, toMask, "-999", "Value = 0")
                maskedBand.save(toMask.name[:-4] + "_Masked.tif")
                # arcpy.Delete_management(bandToMask)
                del toMask
    del maskBand

# delete all _SA bands except *lst and *ndvi
def cleanUp(currentFolder):
    if len(str(currentFolder)) > 0:
        arcpy.workspace = currentFolder
        tempRasterList = arcpy.ListRasters("*","TIF")
        for rasters in tempRasterList:
            if "gamma" in rasters or "delta" in rasters or "lSen" in rasters or "tSen" in rasters or "NDVI_Masked.tif" in rasters or "NDVI.tif" in rasters:
                print "Deleting: ",
                print rasters
                arcpy.Delete_management(rasters)
        del tempRasterList

'''   
#delete any previously-made rasters -- careful, this will delete any of the JMS intermediate rasters, as well as the LST
for root, dirs, files in os.walk(tileFolderRoot):
    try:
        deleteRasters(root, files)
    except:
        print "Failed to delete rasters in: ",
        print str(root)
        error_list = open(error_listLocation, "a")
        error_list.write("File deletion failure in: " + str(root))
        error_list.close()
 '''       

# mask the rasters with the cloudmask


def mask_imgs(img_folder):
    for root, dirs, files in os.walk(img_folder):
        try:
            maskBands(root, files)
        except:
            print "Failed to mask images: ",
            print str(files)
            error_list = open(error_listLocation, "a")
            error_list.write("MaskBand failure in: " + str(files))
            error_list.close()

def extract_sa():
    # Extract the study area raw data based on the 2012 quarantine zone shapefile
    # ACTUALLY I CHANGED THIS TO THE DISSOLVED 5 TOWN SA
    for tile in os.listdir(tileFolderRoot):
        for sensor in os.listdir(os.path.join(tileFolderRoot, tile)):
            for image in os.listdir(os.path.join(tileFolderRoot, tile, sensor)):
                try:
                    env.workspace = os.path.join(tileFolderRoot, tile, sensor, image)
                    imageFolder = os.path.join(tileFolderRoot, tile, sensor, image)
                    bandList = arcpy.ListRasters("*", "TIF")
                    extractedRaster = ""
                    fcTownsSA = "D:\\sync\\gis\\Dissertation\\StudyArea\\StudyAreaDatabase.gdb\\SAExtractionPolygon"
                    # fcQuarZone = "D:\\sync\\gis\\Dissertation\\StudyArea\\StudyAreaDatabase.gdb\\ALBQUarantineWorcester_USDA_20120101"
                    for band in bandList:
                        if len(str(band)) > 0 and "None" not in band:
                            if "_SA" not in band:
                                if "NDVI_Masked" in band or "B6" in band or "B10" in band:
                                    extractedRaster = ExtractByMask(band, fcTownsSA)
                                    extractedRaster.save(str(band)[0:-4] + "_SA.TIF")
                                    del extractedRaster
                except:
                    print "Failed to extract image: ",
                    print str(image)
                    error_list = open(error_listLocation, "a")
                    error_list.write("SA extraction failure in: " + str(image))
                    error_list.close()


def lst():
    # Loop through the sorted imagery folders (tiles) and do all the thermal calculations in sequence
    for tile in os.listdir(tileFolderRoot):
        for sensor in os.listdir(os.path.join(tileFolderRoot, tile)):
            for image in os.listdir(os.path.join(tileFolderRoot, tile, sensor)):
                env.workspace = os.path.join(tileFolderRoot, tile, sensor, image)
                imageFolder = os.path.join(tileFolderRoot, tile, sensor, image)
                os.chdir(imageFolder)
                metadataFile = ""
                thermalBand = ""
                ndviBand = ""
                fileList = os.listdir(imageFolder)
                for file in fileList:
                    if "MTL.txt" in file:
                        metadataFile = file
                    if "LC8" in imageFolder:
                        bandList = arcpy.ListRasters("*_SA*", "TIF")
                        for band in bandList:
                            if "B10_SA" in band and "lSen" not in band and "tSen" not in band and "gamma" not in band and "delta" not in band:
                                thermalBand = arcpy.Raster(band)
                            elif "NDVI_Masked_SA" in band and "emis" not in band :
                                ndviBand = arcpy.Raster(band)
                    elif "LE7" in imageFolder:
                        bandList = arcpy.ListRasters("*_SA*", "TIF")
                        for band in bandList:
                            if "B6_VCID_1_SA" in band and "lSen" not in band and "tSen" not in band and "gamma" not in band and "delta" not in band:
                                if "VCID_1" in band:
                                    thermalBand = arcpy.Raster(band)
                            elif "NDVI_Masked_SA" in band and "emis" not in band:
                                ndviBand = arcpy.Raster(band)
                    elif "LT5" in imageFolder:
                        bandList = arcpy.ListRasters("*_SA*", "TIF")
                        for band in bandList:
                            if "B6_SA" in band and "lSen" not in band and "tSen" not in band and "gamma" not in band and "delta" not in band:
                                thermalBand = arcpy.Raster(band)
                            elif "NDVI_Masked_SA" in band and "emis" not in band:
                                ndviBand = arcpy.Raster(band)
                if len(str(thermalBand)) > 0:
                    if "012030" not in str(thermalBand):
                        try:
                            # print "NDVI BAND IS: ",
                            # print ndviBand
                            lSenRaster = lSen(imageFolder, thermalBand, metadataFile)
                            tSenRaster = tSen(imageFolder, lSenRaster)
                            gammaDeltaList = gammaDelta(imageFolder, tSenRaster, lSenRaster)
                            emisRaster = emissivity(imageFolder, ndviBand)
                            waterVapor = w(thermalBand)
                            psiValues = psi123(radiosoundDict, waterVapor, thermalBand)
                            print "Calculating lst for : ",
                            print thermalBand
                            lstRaster = lst(imageFolder, lSenRaster, gammaDeltaList, emisRaster, psiValues)
                            # cleanUp(imageFolder)
                        except arcpy.ExecuteError:
                            print arcpy.AddError(arcpy.GetMessages(2))
                            print "Could not calculate LST for: ",
                            print imageFolder
                            error_list = open(error_listLocation, "a")
                            error_list.write("LST calculation fail outside module: " + str(imageFolder))
                            error_list.close()


endTime = datetime.datetime.now()
print("Processing completed: ", datetime.datetime.now(), "Total processing time: ", endTime - startTime)