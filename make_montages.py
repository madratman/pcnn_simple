import os

data_dir = '/home/ratneshmadaan/data/turkey_512/Ground_Truth_of_Powerline_Database_512x512/Visible Light (VL)/VL_Original (VL_ORG)'
data_dir_call = '/home/ratneshmadaan/data/turkey_512/Ground_Truth_of_Powerline_Database_512x512/Visible\ Light\ \(VL\)/VL_Original\ \(VL_ORG\)'

for each_image in sorted(os.listdir(data_dir)):
	if 'pcnn' not in each_image and '.DS_Store' not in each_image:
		print each_image
		os.system('montage  ' + os.path.join(data_dir_call, each_image))

montage car*.png -tile 5x1 -geometry +1+1 montage_car.png
montage ultra*.png -tile 5x1 -geometry +1+1 montage_usseq.png
