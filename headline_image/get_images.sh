# Catalogue is too big for github

# TODO replace with wget command when on ADC
freqs="072-103 103-134 139-170 170-231"

for freq in $freqs
do

    if [[ ! -e IDR1_XG_${freq}MHz_rescaled.fits ]]
    then
        scp nhurleywalker@magnus.pawsey.org.au:/astro/mwasci/nhurleywalker/GLEAM-X_headline/headline_${freq}MHz.fits ./
    fi

done

