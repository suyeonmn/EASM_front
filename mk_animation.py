import sys
from PIL import Image, ImageDraw
import glob

# Create the frames
frames = []
#dir = '/data2/utsumi/data/GridFront/v1.0/201006/fig/'
dir = '/home/suyeon/front/copy/'
#imgs = glob.glob(dir+"*.png")
imgs = sorted( glob.glob(dir+"*.png"))
for i in imgs:
        print(i, imgs)
        new_frame = Image.open(i)
        frames.append(new_frame)
#sys.exit()

# Save into a GIF file that loops forever
frames[0].save('./Front_2010_June16-July15.gif', format='GIF', append_images=frames[:], save_all=True, duration=300, Loop=0)
#frames[0].save('./Front_2010_June.gif', format='GIF', append_images=frames[1:], save_all=True, duration=300, Loop=0)




#paths = [ Image.open(i) for i in path]
#imageio.mimsave('./test.gif', paths, fps=0.5)


