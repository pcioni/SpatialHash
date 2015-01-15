#include <cstdlib>
#include <cmath>
#include <vector>
#include "image.h"




// ===================================================================================================
// ===================================================================================================


//function finds the next bucket with the most members in vector mappings
void find_max(std::vector<std::vector<int> >& mappings, int& i_out, int& j_out){
	int max = 0;
	for (int i = 0; i < mappings.size(); ++i) {
		for (int j = 0; j < mappings[i].size(); ++j) {
			if (mappings[i][j] >= max) {
				max = mappings[i][j];
				i_out = i;
				j_out = j;
				
			}
		}
	}
	mappings[i_out][j_out] = 0;
}

//function to make our mappings vector
void get_mappings(const Image<Color> &input, const int& s_offset, std::vector<std::vector<int> > &mappings) {
	for (int i=0;i<input.Width();++i){
		for (int j=0;j<input.Height();++j){
			if (!input.GetPixel(i,j).isWhite())
				++mappings[i%s_offset][j%s_offset];
		}
	}
}

//function to check if we still have non-empty buckets to hash
bool remaining(const std::vector<std::vector<int> > &mappings) {
	for (int i = 0; i < mappings.size(); ++i) {
		for (int j = 0; j < mappings[i].size(); ++j) {
			if (mappings[i][j] != 0)
				return true;
		}
	}
	
	return false;
}

void Compress(const Image<Color> &input, 
              Image<bool> &occupancy, Image<Color> &hash_data, Image<Offset> &offset) {
				  
	occupancy.Allocate(input.Width(), input.Height());
	occupancy.SetAllPixels(false);
	int non_white = 0;
	
	//make our occupancy while counting the non-white pixels
	for (int i = 0; i < input.Width(); ++i){
		for (int j = 0; j < input.Height(); ++j){
			if (!input.GetPixel(i, j).isWhite()) {
				occupancy.SetPixel(i, j, true);
				++non_white;
			}
		}
	}
	
	//functions provided in the HW PDF
	int s_hash = ceil(sqrt(1.01 * non_white));
	int s_offset = ceil(sqrt(non_white / 4));
	
	//make sure s_hash and s_offset don't share a common factor
	if (s_hash%s_offset == 0) 
		s_hash += 1;
	hash_data.Allocate(s_hash, s_hash);
	offset.Allocate(s_offset, s_offset);

	//2d vector of ints to temp store the number of mappings
	std::vector<std::vector<int > > mappings(s_offset,std::vector<int>(s_offset,0));
	get_mappings(input, s_offset, mappings);
	
	//values that will let us retrieve the next most populated offset pixel
	int i_max, j_max; 
	std::vector<std::pair<int,int> > new_mapped;
	
	//loop through all the cells in offset, with most maps to last maps
	while (remaining(mappings)) {
		//set imax and jmax to be the next bucket with the most elements in mappings
		find_max(mappings, i_max, j_max);
		//loop through input pixels
		for (int i = 0; i < input.Width();++i){
			for (int j = 0; j < input.Height(); ++j){
				if (i%s_offset == i_max && j%s_offset == j_max && !input.GetPixel(i,j).isWhite()) {
					//construct a vector of all the pixels that map to the current i_max, j_max
					new_mapped.push_back(std::make_pair(i,j));
				}
			}
		}
		
		int dx_offset = 0;
		int dy_offset = 0;

		while (true) {
			bool good = true;
			for (int m=0; m < new_mapped.size(); ++m) {
				//see if our current dx/dy offsets would cause a collision
				Color tmp = hash_data.GetPixel( (new_mapped[m].first+dx_offset)%s_hash, (new_mapped[m].second+dy_offset )%s_hash ); 
				if ( !tmp.isWhite() )
					good = false;
			}
			if (good) { //if our dx/dy values are good, set our offset and hash values 
				offset.SetPixel(i_max, j_max, Offset(dx_offset, dy_offset));
				for (int i=0; i < input.Width(); ++i) {
					for (int j=0; j < input.Height(); ++j) {
						//assign the hash locations
						if (i%s_offset == i_max && j%s_offset == j_max && !input.GetPixel(i,j).isWhite()) {
							int dx = offset.GetPixel(i%s_offset, j%s_offset).dx;
							int dy = offset.GetPixel(i%s_offset, j%s_offset).dy;
							hash_data.SetPixel((i+dx)%s_hash, (j+dy)%s_hash, input.GetPixel(i,j));
						}
					}
				}
				break;
			}
		
			//base case for if the hash tabel literally won't work.
			//reset our mappings, and make a new hash & offset table
		   else if (dx_offset == 15 && dy_offset == 15) {
				//reset our mappings vector
				for (int i = 0; i < mappings.size(); ++i) {
					for (int j = 0; j < mappings[i].size(); ++j)
						mappings[i][j] = 0;
				}
				get_mappings(input, s_offset, mappings);
				s_hash += 1; //make a new hash table
				hash_data.Allocate(s_hash, s_hash);
				offset.Allocate(s_offset, s_offset);
				dx_offset = 0; //reset dx/dy
				dy_offset = 0;
			}
			//increment our dx and dy one at time, and reset when dx hits the max
			else if (dx_offset == 15) {
				dx_offset = 0;
				++dy_offset;
			} else if (dx_offset < 16) 
				++dx_offset;
				
		} //END OF while(true)
		new_mapped.clear();
		
	} //END OF while(x > 0)
} //END OF FUNCTION


void UnCompress(const Image<bool> &occupancy, const Image<Color> &hash_data, const Image<Offset> &offset, 
                Image<Color> &output) {

	//find offset and hash sizes
	int off_s = offset.Width();  
	int hash_s = hash_data.Width();

	//allocate output size to be equal to occupancy size
	output.Allocate(occupancy.Width(), occupancy.Height());

	//loop through occupancy to find black pixels
	for (int i = 0; i < occupancy.Width(); ++i) {
		for (int j = 0; j < occupancy.Height(); ++j) {
			if (occupancy.GetPixel(i, j)) { //returns true if pixel is black
				Offset tmp = offset.GetPixel(i % off_s, j % off_s); //find dx and dy of the correct offset pixel
				//find the corresponding location in the hash table and set map it to our output
				Color current_hash_loc = hash_data.GetPixel((i + tmp.dx) % hash_s, (j + tmp.dy) % hash_s);
				output.SetPixel(i, j, current_hash_loc); 
			}
		}
	}
}


// ===================================================================================================
// ===================================================================================================

// Takes in two 24-bit color images as input and creates a b&w output
// image (black where images are the same, white where different)
void Compare(const Image<Color> &input1, const Image<Color> &input2, Image<bool> &output) {

  // confirm that the files are the same size
  if (input1.Width() != input2.Width() ||
      input1.Height() != input2.Height()) {
    std::cerr << "Error: can't compare images of different dimensions: " 
         << input1.Width() << "x" << input1.Height() << " vs " 
         << input2.Width() << "x" << input2.Height() << std::endl;
  } else {
    // make sure that the output images is the right size to store the
    // pixel by pixel differences
    output.Allocate(input1.Width(),input1.Height());
    int count = 0;
    for (int i = 0; i < input1.Width(); i++) {
      for (int j = 0; j < input1.Height(); j++) {
        Color c1 = input1.GetPixel(i,j);
        Color c2 = input2.GetPixel(i,j);
        if (c1.r == c2.r && c1.g == c2.g && c1.b == c2.b)
          output.SetPixel(i,j,true);
        else {
          count++;
          output.SetPixel(i,j,false);
        }
      }      
    }     

    // confirm that the files are the same size
    if (count == 0) {
      std::cout << "The images are identical." << std::endl;
    } else {
      std::cout << "The images differ at " << count << " pixel(s)." << std::endl;
    }
  }
}

// ===================================================================================================
// ===================================================================================================

// to allow visualization of the custom offset image format
void ConvertOffsetToColor(const Image<Offset> &input, Image<Color> &output) {
  // prepare the output image to be the same size as the input image
  output.Allocate(input.Width(),input.Height());
  for (int i = 0; i < output.Width(); i++) {
    for (int j = 0; j < output.Height(); j++) {
      // grab the offset value for this pixel in the image
      Offset off = input.GetPixel(i,j);
      // set the pixel in the output image
      int dx = off.dx;
      int dy = off.dy;
      assert (dx >= 0 && dx <= 15);
      assert (dy >= 0 && dy <= 15);
      // to make a pretty image with purple, cyan, blue, & white pixels:
      int red = dx * 16;
      int green = dy *= 16;
      int blue = 255;
      assert (red >= 0 && red <= 255);
      assert (green >= 0 && green <= 255);
      output.SetPixel(i,j,Color(red,green,blue));
    }
  }
}

// ===================================================================================================
// ===================================================================================================

void usage(char* executable) {
  std::cerr << "Usage:  4 options" << std::endl;
  std::cerr << "  1)  " << executable << " compress input.ppm occupancy.pbm data.ppm offset.offset" << std::endl;
  std::cerr << "  2)  " << executable << " uncompress occupancy.pbm data.ppm offset.offset output.ppm" << std::endl;
  std::cerr << "  3)  " << executable << " compare input1.ppm input2.ppm output.pbm" << std::endl;
  std::cerr << "  4)  " << executable << " visualize_offset input.offset output.ppm" << std::endl;
}

// ===================================================================================================
// ===================================================================================================

int main(int argc, char* argv[]) {
  if (argc < 2) { usage(argv[1]); exit(1); }

  if (argv[1] == std::string("compress")) {
    if (argc != 6) { usage(argv[1]); exit(1); }
    // the original image:
    Image<Color> input;
    // 3 files form the compressed representation:
    Image<bool> occupancy;
    Image<Color> hash_data;
    Image<Offset> offset;
    input.Load(argv[2]);
    Compress(input,occupancy,hash_data,offset);
    // save the compressed representation
    occupancy.Save(argv[3]);
    hash_data.Save(argv[4]);
    offset.Save(argv[5]);

  } else if (argv[1] == std::string("uncompress")) {
    if (argc != 6) { usage(argv[1]); exit(1); }
    // the compressed representation:
    Image<bool> occupancy;
    Image<Color> hash_data;
    Image<Offset> offset;
    occupancy.Load(argv[2]);
    hash_data.Load(argv[3]);
    offset.Load(argv[4]);
    // the reconstructed image
    Image<Color> output;
    UnCompress(occupancy,hash_data,offset,output);
    // save the reconstruction
    output.Save(argv[5]);
  
  } else if (argv[1] == std::string("compare")) {
    if (argc != 5) { usage(argv[1]); exit(1); }
    // the original images
    Image<Color> input1;
    Image<Color> input2;    
    input1.Load(argv[2]);
    input2.Load(argv[3]);
    // the difference image
    Image<bool> output;
    Compare(input1,input2,output);
    // save the difference
    output.Save(argv[4]);

  } else if (argv[1] == std::string("visualize_offset")) {
    if (argc != 4) { usage(argv[1]); exit(1); }
    // the 8-bit offset image (custom format)
    Image<Offset> input;
    input.Load(argv[2]);
    // a 24-bit color version of the image
    Image<Color> output;
    ConvertOffsetToColor(input,output);
    output.Save(argv[3]);

  } else {
    usage(argv[0]);
    exit(1);
  }
}

// ===================================================================================================
// ===================================================================================================
