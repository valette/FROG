import argparse
import nibabel as nib

def flip_image( filename, output = "output.nii.gz", orientation = "RAS", threshold = None, fake_orientation = None ):

  image_obj = nib.load( filename )
  if orientation == None : orientation = "RAS"
  current_orientation = nib.aff2axcodes( image_obj.affine )
  print( "Current orientation is : ", current_orientation )

  if not fake_orientation and ( current_orientation == tuple( orientation ) ):
    return filename

  else :
    img_orentation = nib.orientations.io_orientation( image_obj.affine )
    wanted_orientation = nib.orientations.axcodes2ornt( orientation )
    transform = nib.orientations.ornt_transform( img_orentation, wanted_orientation )
    flipped = image_obj.as_reoriented( transform )
    affine = flipped.affine
    if fake_orientation:
      wanted_orientation2 = nib.orientations.axcodes2ornt( fake_orientation )
      transform2 = nib.orientations.ornt_transform( img_orentation, wanted_orientation2 )
      flipped2 = image_obj.as_reoriented( transform2 )
      affine = flipped2.affine
      print( "Real orientation is : ", nib.aff2axcodes( flipped.affine ) )


    #Set Qcode to 1 that the Qform matrix can be used into the further processing
    flipped.header['qform_code'] = 1
    img_data = flipped.get_fdata()
    if threshold : img_data[ img_data < threshold ] = 0
    new_img = img_data.astype( flipped.header.get_data_dtype() )
    img_conv = nib.Nifti1Image( new_img, affine, flipped.header)
    nib.save( img_conv, output )      
    print("New orientation is : ", nib.aff2axcodes( img_conv.affine ))
    return output


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Flip NIFTI image')
    parser.add_argument('--orientation', help = "Output volume orientation" )
    parser.add_argument('--fake_orientation', help = "Fake output volume orientation" )
    parser.add_argument('--threshold', help = "Threshold", type = float )
    parser.add_argument('input_volume', help = "Input volume" )
    parser.add_argument('--output', help = "Output file name", default = "output.nii.gz" )
    args = parser.parse_args()
    flip_image( args.input_volume, output = args.output, orientation = args.orientation, fake_orientation = args.fake_orientation, threshold = args.threshold )
