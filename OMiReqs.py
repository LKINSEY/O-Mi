#%% IMPORTS
## IMPORTS
import requests, sys, json, random, os,  s3fs, zarr
import multiprocessing as mp
from pathlib import Path
from cloudvolume import CloudVolume
from ftplib import FTP
import numpy as np
import hyperspy.api as hs
from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.by import By
import tifffile as tf

def readData(url, urlDir='C:/temp',crop=(500,500,500)):

    saveDir = Path(urlDir) / 'data'
    saveDir.mkdir(parents=True, exist_ok=True)

    if 'idr.openmicroscopy.org' in url: #size input
        from omero.gateway import BlitzGateway
        methodUsed = 'omero-py'
        conn = BlitzGateway('public', 'public', host='idr.openmicroscopy.org', secure=True)
        conn.connect()
        print('Open Microscope OMERO API URL:\n', url)
        imID = url[url.find('img_detail/')+11:url.find('/?dataset')]
        IM = conn.getObject("Image", imID)
        metadata = requests.get(f'https://idr.openmicroscopy.org/webclient/imgData/{imID}/').json()
        print('Grabbing:\n', IM.getName(), f'\nDimensions: {IM.getSizeZ()}x{IM.getSizeY()}x{IM.getSizeX()}x{IM.getSizeC()}x{IM.getSizeT()}')
        if crop[0]*2 <IM.getSizeZ():
            cropping = True
            z = random.randint(crop[0]+1,IM.getSizeZ()-crop[0]-1)
            zBounds = (z-crop[0], z+crop[0])
            print('NOTICE: Cropping z:', crop[0]*2)
        else:
            zBounds = (0, IM.getSizeZ())
        if crop[1]*2 < IM.getSizeY():
            cropping = True
            y = random.randint(crop[1]+1,IM.getSizeY()-crop[1]-1)
            yBounds = (y-crop[1], y+crop[1])
            print('NOTICE: Cropping y:', crop[1]*2)
        else:
            yBounds = (0, IM.getSizeY())
        if crop[2]*2 < IM.getSizeX():
            cropping = True
            x = random.randint(crop[2]+1,IM.getSizeX()-crop[2]-1)
            xBounds = (x-crop[2], x+crop[2])
            print('NOTICE: Cropping x:', crop[2]*2)
        else:
            xBounds = (0, IM.getSizeX())
        if cropping:
            print('NOTICE: Cropping size to shape:', (zBounds[1]-zBounds[0], yBounds[1]-yBounds[0], xBounds[1]-xBounds[0]))
        else: 
            print(f'NOTICE: Crop size {crop} is too large for image dimensions {IM.getSizeZ()}x{IM.getSizeY()}x{IM.getSizeX()}. No Crop Necessary.')
        vol = np.ndarray([zBounds[1]-zBounds[0], yBounds[1]-yBounds[0], xBounds[1]-xBounds[0]], dtype=np.uint16) 
        count = 0
        for i in range(zBounds[0], zBounds[1]):
            print(f'\rGrabbing plane {i+1} of {IM.getSizeZ()}', end='', flush=True)
            plane = IM.getPrimaryPixels().getPlane(i,0,0) 
            vol[count,:,:] = plane[yBounds[0]:yBounds[1], xBounds[0]:xBounds[1]]
            count += 1

        meta={
            'Name': IM.getName(),
            'imID': imID,
            'SizeX': xBounds[1]-xBounds[0],
            'SizeY': yBounds[1]-yBounds[0],
            'SizeZ': zBounds[1]-zBounds[0],
            'SizeC': IM.getSizeC(),
            'SizeT': IM.getSizeT(),
            'method': methodUsed,
            'url': url,
            'info': metadata
        }
        conn.close()

    elif 'neuroglancer' in url: #size input but default 1000x1000x1000 
        methodUsed = 'cloud-volume'
        print('Hemibrain Cloud-Volume API:\n', url)
        # block = 500
        source = f'https://storage.googleapis.com/'+url[url.find('precomputed://gs://neuroglancer')+19:url.find('/jpeg')+5]
        print('Grabbing:\n',source)
        data = CloudVolume(
            source,
            mip=0, 
            use_https=True,
            parallel=True,
            secrets=None
        )

        z,y,x = data.info['scales'][0]['size']
        z,y,x = [random.randint(crop[0]+1, z-crop[0]-1), random.randint(crop[1]+1, y-crop[1]-1), random.randint(crop[2]+1, x-crop[2]-1)]

        vol = data[z-crop[0]:z+crop[0], y-crop[1]:y+crop[1], x-crop[2]:x+crop[2]] 

        meta={
            'Name': source,
            'imID': (url[url.find('precomputed://gs://neuroglancer')+19:url.find('/jpeg')]).replace('/','_'),
            'SizeX': x,
            'SizeY': y,
            'SizeZ': z,
            'SizeC': None,
            'SizeT': None,
            'method': methodUsed,
            'url': url,
            'info': data.info
        }

    elif 'ebi' in url: #no size input
        print('EMPIAR API:\n', url)
        methodUsed = 'FTP'
        entryID = url[url.find('EMPIAR-')+7:-1]
        source = 'https://www.ebi.ac.uk/empiar/api/entry/'+entryID
        print('Entry ID:\n', entryID)
        metadata = requests.get(source).json()
        z = int(metadata['EMPIAR-11759']['imagesets'][0]['num_images_or_tilt_series'])
        x = int(metadata['EMPIAR-11759']['imagesets'][0]['image_width'])
        y = int(metadata['EMPIAR-11759']['imagesets'][0]['image_height'])
        
        vol = np.ndarray([z,x,y], dtype=np.uint8) 

        #access open source database via FTP and grab data
        ftp = FTP('ftp.ebi.ac.uk')
        ftp.login()
        ftp.cwd(f"/empiar/world_availability/{entryID}/data")
        files = ftp.nlst()

        #for now hardcoded to only read DM3 files
        for i, f in enumerate(files):
            print(f)
            if f.endswith(metadata['EMPIAR-11759']['imagesets'][0]['data_format'].lower()):
                t=i-1 #first file is xml metadata
                with open(f, 'wb') as localFile:
                    ftp.retrbinary('RETR ' + f, localFile.write)
                # attempting to read DM3 file with hyperspy, saves temp local im then deletes it after writing to arr
                im = hs.load(f)
                arr = im.data
                ax, ay = arr.shape
                if ax < x or ay < y:
                    print(f'WARNING: image dimensions {ax}x{ay} do not match expected dimensions {x}x{y}')
                    newX = x-ax
                    newY = y-ay
                    arr = np.pad(arr, ((0, newX), (0, newY)), mode='constant', constant_values=np.nan)
                vol[t,:,:] = arr
                os.remove(f)

        vol = vol[:,:ax, :ay] #crop to actual dimensions if padding was added
        if crop[0]*2 < z:
            cropping = True
            z = random.randint(crop[0]+1, z-crop[0]-1)
            zBounds = (z-crop[0], z+crop[0])
        else:
            zBounds = (0, z)
        if crop[1]*2 < ay:
            cropping = True
            y = random.randint(crop[1]+1, ay-crop[1]-1)
            yBounds = (y-crop[1], y+crop[1])
        else:
            yBounds = (0, ay)
        if crop[2]*2 < ax:
            cropping = True
            x = random.randint(crop[2]+1, ax-crop[2]-1)
            xBounds = (x-crop[2], x+crop[2])
        else:
            xBounds = (0, ax)
            print('NOTICE: Cropping size to shape:', crop)
        if cropping:
            print('NOTICE: Cropping size to shape:', [zBounds[1]-zBounds[0], yBounds[1]-yBounds[0], xBounds[1]-xBounds[0]])
            vol = vol[zBounds[0]:zBounds[1], yBounds[0]:yBounds[1], xBounds[0]:xBounds[1]]

        meta={
            'Name': source,
            'imID': 'EMPIAR-' + entryID,
            'SizeX': ax,
            'SizeY': ay,
            'SizeZ': z,
            'SizeC': None,
            'SizeT': None,
            'method': methodUsed,
            'url': url,
            'info': (metadata)
        }
        ftp.quit()

    elif 'epfl' in url: #no size input
        print('EPFL URL:\n', url)
        print('NOTICE: No size information available, downloading entire dataset from download page. This may take a while...')
        print('Pretending to run this because I want to save time...')
        methodUsed = 'Webpage'  
        source = url
        imID = (url[url.find('/cvlab/'):]).replace('/','_')

        # #avoid website's anti-scraping security using simulated chrome
        # driver = webdriver.Chrome()
        # driver.get(url)
        # elem = driver.find_element(By.XPATH, '//a[contains(text(),"Download dataset")]')
        # download_url = elem.get_attribute("href")
        # print('Download URL:\n', download_url)
        # driver.quit()

        # #data download
        # with requests.get(download_url, stream=True) as r:
        #     r.raise_for_status()
        #     total = int(r.headers.get('content-length', 0))
        #     with open(urlDir+f'data/{imID}/{imID}_data.zip', 'wb') as f, tqdm(
        #         total = total, unit='B', unit_scale=True, desc=download_url.split('/')[-1]
        #     ) as pbar:
        #         for chunk in r.iter_content(chunk_size=8192):
        #             if chunk:
        #                 f.write(chunk)
        #                 pbar.update(len(chunk))
        

        meta={
            'Name': source,
            'imID': imID,
            'SizeX': None,
            'SizeY': None,
            'SizeZ': None,
            'SizeC': None,
            'SizeT': None,
            'method': methodUsed,
            'url': url,
            'info': None
        } 
        vol = None

    elif 'openorganelle' in url: #size input
        #TODO: find way to differentiate zarr from N5 metadata, currently hardcoded to read zarr
        print('OpenOrganelle URL:\n', url)
        methodUsed = 'S3-zarr'
        imID = url[url.find('datasets/')+9:]
        source = f'janelia-cosem-datasets/{imID}/{imID}.zarr'
        fs = s3fs.S3FileSystem(anon=True)
        group = zarr.open(s3fs.S3Map(root=source,s3=fs, check=False), mode='r')
        
        
        data = group['recon-2/em/fibsem-int16/s0']
        #must call data to get chunk
        z, y, x = data.shape
        if crop[0]*2 < z:
            cropping = True
            z = random.randint(crop[0]+1, z-crop[0]-1)
            zBounds = (z-crop[0], z+crop[0])
            print('NOTICE: Cropping z:', crop[0]*2)
        else:
            zBounds = (0, z)
        if crop[1]*2 < y:
            cropping = True
            y = random.randint(crop[1]+1, y-crop[1]-1)
            yBounds = (y-crop[1], y+crop[1])
            print('NOTICE: Cropping y:', crop[1]*2)
        else:
            yBounds = (0, y)
        if crop[2]*2 < x:
            cropping = True
            x = random.randint(crop[2]+1, x-crop[2]-1)
            xBounds = (x-crop[2], x+crop[2])
            print('NOTICE: Cropping x:', crop[2]*2)
        else:
            xBounds = (0, x)

        if cropping:
            print('NOTICE: Cropping size to shape:', [zBounds[1]-zBounds[0], yBounds[1]-yBounds[0], xBounds[1]-xBounds[0]])
        vol = data[zBounds[0]:zBounds[1], yBounds[0]:yBounds[1], xBounds[0]:xBounds[1]]

        # print(dir(group))
        meta={
            'Name': source,
            'imID': imID,
            'SizeX': data.shape[2],
            'SizeY': data.shape[1],
            'SizeZ': data.shape[0],
            'SizeC': None,
            'SizeT': None,
            'method': methodUsed,
            'url': url,
            'info': str(group.info)
        }
    
    dataIDDIR = Path(saveDir) / f'{meta["imID"]}'
    dataIDDIR.mkdir(parents=True, exist_ok=True)
    with open(str(dataIDDIR)+'/metadata.json', 'w') as f:
        json.dump(meta, f, indent=4)
    if vol is not None:
        tf.imwrite(str(dataIDDIR)+f'/{meta["imID"]}.tif', vol)

# '''
# Example stand-alone usage:
# readData(
#     'https://openorganelle.janelia.org/datasets/jrc_mus-nacc-2',
#     urlDir='D:/openMicro',
#     crop=(10,256,256)
# )
# '''

#%% IDENTIFY URL
## IDENTIFY URL

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Usage: python OMiReqs.py <path_to_URL_directory>')
        sys.exit(1)
    else:
        # print(sys.argv[1])
        urlDir = sys.argv[1] #path to directory with .txt files containing URLs specified in CLI
        print('Searching for Stored URLs in:\n', urlDir)
        
        with open(urlDir+'/imageURLs.txt') as f:
            urls = f.read().splitlines()
        print('URLs:\n', urls)

        args = [(url, urlDir, (10,256,256)) for url in urls] #can add CLI inputs here in future
        with mp.Pool(processes=mp.cpu_count()) as pool:
            metaData = pool.starmap(readData, args)

    
