import os

data_dir = '/home/ratneshmadaan/data/turkey_512/Ground_Truth_of_Powerline_Database_512x512/Visible Light (VL)/VL_Original (VL_ORG)'
data_dir_call = '/home/ratneshmadaan/data/turkey_512/Ground_Truth_of_Powerline_Database_512x512/Visible\ Light\ \(VL\)/VL_Original\ \(VL_ORG\)'

for each_image in sorted(os.listdir(data_dir)):
	if 'pcnn' in each_image and '.DS_Store' not in each_image:
		print each_image
		img_orig = cv2.imread(os.path.join(data_dir, each_image))
		img_basename = img_path.split('/')[-2].split('.')[0] 

		# todo some morphology. try cv2.getcontours ? (make opencv 3 docker image)

		lines = cv2.HoughLines(img_canny, 1, np.pi/180, 100)

		if lines is not None:
			for rho,theta in lines[0]:
				a = np.cos(theta)
				b = np.sin(theta)
				x0 = a*rho
				y0 = b*rho
				x1 = int(x0 + 1000*(-b))
				y1 = int(y0 + 1000*(a))
				x2 = int(x0 - 1000*(-b))
				y2 = int(y0 - 1000*(a))
				cv2.line(image_clone,(x1,y1),(x2,y2),(0,0,255),2)
			
		# cv2.imwrite('results/' + img_basename + '_' + str(idx) + '.png', np.concatenate([each_window, image_clone], axis=1))

		# cv2.imwrite('results/' + img_basename + '_' + str(idx) + '.png', np.concatenate([img_orig, img_canny_montage, img_lines], axis=1))