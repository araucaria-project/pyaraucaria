from astropy.io import fits
import numpy


def save_fits_from_array(array, file_name, header):

    file_name = f"{file_name}.fits"

    hdr = fits.Header()
    if isinstance(header, dict):
        for n in header.key():
            hdr[n] = header[n]
    else:
        hdr["FITS_STD"] = "No fits header loaded"

    hdu = fits.PrimaryHDU(data=array, header=hdr)
    hdul = fits.HDUList([hdu])
    hdul.writeto(file_name)



"""
    hdr["FITS_STD"]="BETA1"
    hdr["TEL"]=str(self.parent.active_tel)
    hdr["UT"]=str(self.parent.ut)
    hdr["JD"]=str(self.parent.jd)

    hdr["TEL_RA"]=str(self.parent.mount_ra )
    hdr["TEL_DEC"]=str(self.parent.mount_dec)
    hdr["TEL_ALT"]=str(self.parent.mount_alt)
    hdr["TEL_AZ"]=str(self.parent.mount_az )

    hdr["TYPE"]=self.parent.ob_type
    hdr["OBJECT"]=self.parent.ob_name
    hdr["FILTER"]=self.parent.curent_filter
    hdr["DIT"]=self.parent.dit_exp

    hdr["CCD_NAME"]=str(self.parent.ccd_name)
    hdr["CCD_TEMP"]=str(self.parent.ccd_temp)
    hdr["CCD_BINX"]=str(self.parent.ccd_binx)
    hdr["CCD_BINY"]=str(self.parent.ccd_biny)
    hdr["RAED_MOD"]=str(self.parent.ccd_readmode)


    hdr["FOCUS"]=str(self.parent.focus_value)"""