import numpy as np
import pylab as pyl
from matplotlib.pyplot import ion
import pyfits as pf
from astropy.io import fits
from scipy import *
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import Distance, SkyCoord


ion()

###
# RUN USING
# aegean \
# --seedclip=4 \
# --maxsummits=5 \
# --cores 1 \
# --autoload \
# --progress \
# --nopositive \
# --negative \
# --region /astro/mwasci/tvernstrom/GLEAMX/negative/prime_plus_edge.mim \
# --psf=/astro/mwasci/tgalvin/Mosaic_GLEAMX_Version2/Coadd_Mosaic/IDR_v1/DeepAegeanFix/XG_170-231MHz_projpsf_psf.fits \
# --table=/astro/mwasci/tvernstrom/GLEAMX/negative/XG_170-231MHz_psf_comp_negative_dafix1.fits \
# /astro/mwasci/tgalvin/Mosaic_GLEAMX_Version2/Coadd_Mosaic/IDR_v1/DeepAegeanFix/XG_170-231MHz.fits

cat_neg=fits.getdata('XG_170-231MHz_psf_comp_negative_dafix_comp.fits')

cat_full=fits.getdata('XG_170-231MHz_comp_rescaled.fits',1)

hd=fits.getheader('XG_170-231MHz.fits')


nsamps=((((13.-4.)*15*np.cos(np.radians(-26.7)))*12)/(hd['BMAJ']*hd['BMIN']*1.13))*2.

## FIND THOSE WITHIN THE CENTRAL REGION
##USING MIN RA OF 60 AND MAX 195 AND MIN DEC -32.7 MAX DEC -20.7
keepneg=np.argwhere((cat_neg['ra']>=60)&(cat_neg['ra']<=195)&(cat_neg['dec']>=-32.7)&(cat_neg['dec']<=-20.7))[:,0]

keeppos=np.argwhere((cat_full['ra']>=60)&(cat_full['ra']<=195)&(cat_full['dec']>=-32.7)&(cat_full['dec']<=-20.7))[:,0]


cat_neg=cat_neg[keepneg]
cat_full=cat_full[keeppos]


#CUT OUT ONES BELOW SNR OF 5
cat_neg=cat_neg[abs(cat_neg['int_flux']/cat_neg['local_rms'])>=5]
cat_full=cat_full[abs(cat_full['int_flux']/cat_full['local_rms'])>=5]


nneg=cat_neg.shape[0]
npos=cat_full.shape[0]

noise1=np.median( cat_full['local_rms'])



## MAKE A REGION FILE FOR INSPECTION ##
outf1=open('negatives.reg','w')
print>>outf1,'fk5'

for i in range(cat_neg.size):
    print >> outf1,'point('+str(cat_neg['ra'][i])+','+str(cat_neg['dec'][i])+') # point=cross color=blue'

outf1.close()


outf1=open('positives.reg','w')
print>>outf1,'fk5'

for i in range(cat_full.size):
    print >> outf1,'point('+str(cat_full['ra'][i])+','+str(cat_full['dec'][i])+') # point=X color=red'

outf1.close()


## LOOK FOR NEGATIVES SOURCES NEAR POSTIVE ONES
coords_neg = SkyCoord(cat_neg['ra']*u.deg, cat_neg['dec']*u.deg, frame='fk5')
coords_pos=SkyCoord(cat_full['ra']*u.deg,cat_full['dec']*u.deg,frame='fk5')
search=coords_neg.search_around_sky(coords_pos,15*u.arcmin)


## FIND THE BRIGHTEST POSTIVE SOURCE WITHIN 15 ARCMINUTE OF NEGATIVE AND THE CLOSEST POSTIVE
brightid=np.zeros(nneg)
sep_bright=np.zeros(nneg)


closeid=np.zeros(nneg)
sep_close=np.zeros(nneg)
for i in range(nneg):
    jj=np.argwhere(search[1]==i)[:,0]
    brightid[i]=search[0][jj[np.argmax(cat_full['peak_flux'][search[0][jj]])]]
    sep_bright[i]=search[2][jj[np.argmax(cat_full['peak_flux'][search[0][jj]])]].value*60.
    sep_close[i]=(search[2][jj].value*60).min()
    closeid[i]=search[0][jj][np.argmin(search[2][jj])]

brightid=brightid.astype(int)
closeid=closeid.astype(int)
peak_flux_bright=cat_full['peak_flux'][brightid]
drange=abs(peak_flux_bright/cat_neg['peak_flux'])

## SET CRITERIA OF WITHIN 10 ARCMIN, PEAK FLUX OF POSITIVE >= 0.1 JY AND ABS(FLUX RATION) >= 30
#criteria1=np.argwhere((sep_bright<=10)&(drange>=30)&(peak_flux_bright>=0.1))[:,0]
#criteria2=np.argwhere(sep_close<=2)[:,0]

criteria1=np.argwhere((sep_bright<=13)&(drange>=250)&(peak_flux_bright>=1))[:,0]

flag=np.zeros(nneg)
flag[criteria1]=1
#flag[criteria2]=1

good=np.argwhere(flag==0)[:,0]

## DO SAME MATCHING AND FILTERING FOR POSTIVE CATALOGUE
search_pos=coords_pos.search_around_sky(coords_pos,15*u.arcmin)

brightid_pos=np.zeros(npos)
sep_bright_pos=np.zeros(npos)


closeid_pos=np.zeros(npos)
sep_close_pos=np.zeros(npos)
for i in range(npos):
    jj=np.argwhere(search_pos[1]==i)[:,0]
    kj=np.argwhere(search_pos[0][jj]!=i)[:,0]
    jj=jj[kj]
    if len(jj) !=0:
        brightid_pos[i]=search_pos[0][jj[np.argmax(cat_full['peak_flux'][search_pos[0][jj]])]]
        sep_bright_pos[i]=search_pos[2][jj[np.argmax(cat_full['peak_flux'][search_pos[0][jj]])]].value*60.
        sep_close_pos[i]=(search_pos[2][jj].value*60).min()
        closeid_pos[i]=search_pos[0][jj][np.argmin(search_pos[2][jj])]

brightid_pos=brightid_pos.astype(int)
closeid_pos=closeid_pos.astype(int)
peak_flux_bright_pos=cat_full['peak_flux'][brightid_pos]
drange_pos=abs(peak_flux_bright_pos/cat_full['peak_flux'])



#criteria1_pos=np.argwhere((sep_bright_pos<=10)&(drange_pos>=350)&(peak_flux_bright_pos>=2.))[:,0]

criteria1_pos=np.argwhere((sep_bright_pos<=4)&(drange_pos>=350)&(peak_flux_bright_pos>=2.))[:,0]
criteria2_pos=np.argwhere((sep_bright_pos<=12)&(sep_bright_pos>5)&(drange_pos>=650)&(peak_flux_bright_pos>=6.))[:,0]


flag_pos=np.zeros(npos)
flag_pos[criteria1_pos]=1
flag_pos[criteria2_pos]=1


good_pos=np.argwhere(flag_pos==0)[:,0]


## MAKE A NEW REGION FILE 
outf1=open('negatives_sub_filtered_new.reg','w')
print>>outf1,'fk5'


for i in range(good.size):
    print >> outf1,'point('+str(cat_neg['ra'][good[i]])+','+str(cat_neg['dec'][good[i]])+') # point=X color=green'

outf1.close()

outf1=open('positives_sub_filtered_new.reg','w')
print>>outf1,'fk5'


for i in range(good_pos.size):
    print >> outf1,'point('+str(cat_full['ra'][good_pos[i]])+','+str(cat_full['dec'][good_pos[i]])+') # point=cross color=green'

outf1.close()

cat_neghd=fits.getheader('XG_170-231MHz_psf_comp_negative_dafix_comp.fits',1)
cat_fullhd=fits.getheader('XG_170-231MHz_comp_rescaled.fits',1)

fits.writeto('XG_170-231MHz_comp_rescaled_filtered.fits',cat_full[good_pos],header=cat_fullhd)

fits.writeto('XG_170-231MHz_psf_comp_negative_dafix_comp_filtered.fits',cat_neg[good],header=cat_neghd)


outf1=open('Filtered_Cat_GOOD_IDS.dat','w')
for i in range(good_pos.size):
    print >> outf1,good_pos[i]

outf1.close()

###CALCULATE GAUSSIAN DIST AND MAKE SOME PLOTS

x_full= np.linspace(-8, 8, 1501)
x_full1= (x_full[1:]+x_full[:-1])/2.
area=x_full[1]-x_full[0]
gauss1=exp((-1.)*(x_full1**2)/(2*(1.**2)))
gauss11=(gauss1*area/(np.pi*2)**.5)


fig1=pyl.figure(1,figsize=(9,9))
ax1=pyl.subplot(111)
ax1.tick_params(labelsize=13)
hxy=pyl.hist(cat_neg['int_flux']/cat_neg['local_rms'],bins=30,alpha=.5,color='blue',label='No Filtering')
hxy2=pyl.hist(cat_neg['int_flux'][good]/cat_neg['local_rms'][good],bins=hxy[1],alpha=0.8,color='r',label='Filtered')

#pyl.plot(x_full1,(gauss11)*(nsamps),'k-',label=r'Noise Gaussian' )
pyl.axis([-10,-4,0.,70.])
pyl.legend(loc=2,prop={'size':13})

pyl.ylabel('No.of Detections',fontsize=16)
pyl.xlabel('Signal to Noise Ratio',fontsize=16)
fig1.tight_layout()
fig1.savefig('GLEAMX_neg.png')


hsnrfull=np.histogram(cat_full['peak_flux'][good_pos]/cat_full['local_rms'][good_pos],bins=np.linspace(5,8,10))
hsnrneg=np.histogram((cat_neg['peak_flux'][good]*-1)/cat_neg['local_rms'][good],bins=np.linspace(5,8,10))
hsnrhx=(hsnrneg[1][1:]+hsnrneg[1][:-1])/2

fig1=pyl.figure(1,figsize=(9,9))
ax1=pyl.subplot(111)
ax1.tick_params(labelsize=13)
pyl.plot(hsnrhx,(hsnrneg[0]/hsnrfull[0].astype('float64'))*100,'r-')
pyl.ylabel('% of Fake Sources',fontsize=16)
pyl.xlabel('ABS(Peak SNR)',fontsize=16)
fig1.tight_layout()
fig1.savefig('GLEAMX_neg-ratio.png')


hsnrfull=np.histogram(cat_full['int_flux'][good_pos]/cat_full['local_rms'][good_pos],bins=np.linspace(5,8,10))
hsnrneg=np.histogram((cat_neg['int_flux'][good]*-1)/cat_neg['local_rms'][good],bins=np.linspace(5,8,10))
hsnrhx=(hsnrneg[1][1:]+hsnrneg[1][:-1])/2

fig1=pyl.figure(1,figsize=(9,9))
ax1=pyl.subplot(111)
ax1.tick_params(labelsize=13)
pyl.plot(hsnrhx,(hsnrneg[0]/hsnrfull[0].astype('float64'))*100,'r-')
pyl.ylabel('% of Fake Sources',fontsize=16)
pyl.xlabel('ABS(Total SNR)',fontsize=16)
fig1.tight_layout()
fig1.savefig('GLEAMX_neg-ratio_int.png')
