import java.awt.image.BufferedImage;
import java.awt.image.ConvolveOp;
import java.awt.image.Kernel;
import java.awt.image.RenderedImage;
import java.awt.Graphics;
import java.awt.Image;
import java.awt.Panel;
import java.awt.image.BufferedImageOp;
import java.awt.image.ColorModel;
import java.io.BufferedReader;
import java.io.File;
import java.awt.image.DataBufferByte;
import java.io.IOException;
import java.io.InputStreamReader;
import javax.imageio.ImageIO;
import javax.swing.JFrame;

public class GaborFilter {

	private double[] orientations;
	private double waveLength;
	private double offset;
	private double aspectRatio;
	private double bandwidth;
	private int width;
	private int height;
	private double[][] kernel;

	public GaborFilter(double waveLength, double[] orientations, double offset, double aspectRatio, double bandwidth, int width, int height) {
		this.aspectRatio = aspectRatio;
		this.bandwidth = bandwidth;
		this.width = width;
		this.height = height;
		this.waveLength = waveLength;
		this.orientations = orientations;
		this.offset = offset;
		this.kernel = new double[width][height];
	}
	
	public double[][] getKernelArray() {
		return this.kernel;
	}

	/**
	* Calculates the Sigma for the given Wave Length and Bandwidth
	*
	* @param waveLength - Wave Length
	* @param bandwidth - Bandwidth
	* @return - Sigma (Deviation)
	*/
	private static double calculateSigma(double waveLength, double bandwidth) {
		return waveLength*Math.sqrt(Math.log(2)/2)*(Math.pow(2, bandwidth) + 1)/((Math.pow(2, bandwidth) - 1)*Math.PI);
	}


	private static double gaborFunction(double x, double y, double sigma, double aspectRatio, double waveLength, double offset) {
		double gaborReal = Math.exp(-(Math.pow(x/sigma, 2) + Math.pow(y*aspectRatio/sigma, 2))/2)*Math.cos(2*Math.PI*x/waveLength + offset);
		//gaborReal = 0;
		double gaborImage = Math.exp(-(Math.pow(x/sigma, 2) + Math.pow(y*aspectRatio/sigma, 2))/2)*Math.sin(2*Math.PI*x/waveLength + offset);
		//gaborImage = 0;

		//double gaborComplex = Math.exp(-(Math.pow(x/sigma, 2) + Math.pow(y*aspectRatio/sigma, 2))/2)*Math.exp(Math.Complex(0,2*Math.PI*x/waveLength + offset));

		return Math.sqrt(Math.pow(gaborReal, 2) + Math.pow(gaborImage, 2));
	}


	/**
	* Returns the ConvolveOp for the Gabor Filter
	*
	* @return - ConvolveOp
	*/
	public ConvolveOp getConvolveOp() {
		return new ConvolveOp(getKernel(), ConvolveOp.EDGE_NO_OP, null);
	}


	/* Implementation of Convolve function
	* 
	* 
	*/
	public static double[][] convertTo2DWithoutUsingGetRGB(BufferedImage image) {

		final byte[] pixels = ((DataBufferByte)image.getRaster().getDataBuffer()).getData();
		final int width = image.getWidth();
		final int height = image.getHeight();
		final boolean hasAlphaChannel = image.getAlphaRaster() != null;

		double[][] result = new double[height][width];
		if (hasAlphaChannel) {
			final int pixelLength = 4;
			for (int pixel = 0, row = 0, col = 0; pixel < pixels.length; pixel += pixelLength) {
				int argb = 0;
				argb += (((int) pixels[pixel] & 0xff) << 24); // alpha
				argb += ((int) pixels[pixel + 1] & 0xff); // blue
				argb += (((int) pixels[pixel + 2] & 0xff) << 8); // green
				argb += (((int) pixels[pixel + 3] & 0xff) << 16); // red
				result[row][col] = argb;
				col++;
				if (col == width) {
					col = 0;
					row++;
				}
			}
		} else {
			final int pixelLength = 3;
			for (int pixel = 0, row = 0, col = 0; pixel < pixels.length; pixel += pixelLength) {
				int argb = 0;
				argb += -16777216; // 255 alpha
				argb += ((int) pixels[pixel] & 0xff); // blue
				argb += (((int) pixels[pixel + 1] & 0xff) << 8); // green
				argb += (((int) pixels[pixel + 2] & 0xff) << 16); // red
				result[row][col] = argb;
				col++;
				if (col == width) {
					col = 0;
					row++;
				}
			}
		}

		return result;
	}


	public static double singlePixelConvolution(double [][] input, int x, int y, double [][] k, int kernelWidth, int kernelHeight){
		double output = 0;
		for(int i=0;i<kernelWidth;++i){
			for(int j=0;j<kernelHeight;++j){
				output = output + (input[x+i][y+j] * k[i][j]);
			}
		}
		return output;
	}

	public static int applyConvolution(int [][] input, int x, int y, double [][] k, int kernelWidth, int kernelHeight){
		int output = 0;
		for(int i=0;i<kernelWidth;++i){
			for(int j=0;j<kernelHeight;++j){
				output = output + (int) Math.round(input[x+i][y+j] * k[i][j]);
			}
		}
		return output;
	}


	/**
	* Takes a 2D array of grey-levels and a kernel and applies the convolution
	* over the area of the image specified by width and height.
	*
	* @param input the 2D double array representing the image
	* @param width the width of the image
	* @param height the height of the image
	* @param kernel the 2D array representing the kernel
	* @param kernelWidth the width of the kernel
	* @param kernelHeight the height of the kernel
	* @return the 2D array representing the new image
	*/
	public static double [][] convolution2D(double [][] input, int width, int height, double [][] kernel, int kernelWidth, int kernelHeight){
		int smallWidth = width - kernelWidth + 1;
		int smallHeight = height - kernelHeight + 1;
		double [][] output = new double [smallWidth][smallHeight];
		for(int i=0;i<smallWidth;++i){
			for(int j=0;j<smallHeight;++j){
				output[i][j]=0;
			}
		}
		for(int i=0;i<smallWidth;++i){
			for(int j=0;j<smallHeight;++j){
				output[i][j] = singlePixelConvolution(input,i,j,kernel,
				kernelWidth,kernelHeight);
				//if (i==32- kernelWidth + 1 && j==100- kernelHeight + 1) System.out.println("Convolve2D: "+output[i][j]);
			}
		}
		return output;
	}


	/**
	* Returns the Kernel for the Gabor filter
	*
	* @return - Kernel
	*/
	public void printKernel(double[][] kernel2) {
		for (int i =0; i<kernel2.length; i++) {
			for (int j=0; j< kernel2[i].length; j++) {
				System.out.print(kernel2[i][j] + " ");
			}
			
		}
	}
	

	public Kernel getKernel() {
		double sigma = calculateSigma(waveLength, bandwidth);
		float[] data = new float[width*height];
		for(int k = 0, x = -width/2; x <= width/2; x++) {
			for(int y = -height/2; y <= height/2; y++) {
				for(double orientation : orientations) {
					double x1 = x*Math.cos(orientation) + y*Math.sin(orientation);
					double y1 = -x*Math.sin(orientation) + y*Math.cos(orientation);
					data[k] += (float)(gaborFunction(x1, y1, sigma, aspectRatio, waveLength, offset));
					
					kernel[x+(width/2)][y+ (height/2)] = (double)data[k];
					//printKernel(kernel);

				}
				k++;
			}
		}
		float sum = 0f;
		for(int i = 0; i < width; i++) {
			for(int j = 0; j < height; j++) {
				sum += data[i*j + j];
			}
		}
		sum /= width*height;
		for(int i = 0; i < width; i++) {
			for(int j = 0; j < height; j++) {
				data[i*j + j] -= sum;
			}
		}
		//printKernel(data);
		return new Kernel(width, height, data);
	}


	/**
	* Filters the bufferedImage using the Gabor filter. If the bufferedImageDestination is not null
	* the bufferedImage is used as the destination
	*
	* @param bufferedImage - buffered image to be used as the source
	* @param bufferedImageDestination - buffered image to be used as the destination
	* @return - the rendered image
	*/
	public RenderedImage filter(BufferedImage bufferedImage, BufferedImage bufferedImageDestination) {
		return getConvolveOp().filter(bufferedImage, bufferedImageDestination);
	}


	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		/*
		JFrame frame=new JFrame("display image");
		Panel panel = new GaborFilter();
		frame.getContentPane().add(panel);
		frame.setSize(500,500);
		frame.setVisible(true);
		*/

		File file = new File("./src/Images/tank3-filtered.jpg");
		File file2 = new File("./src/Images/test.jpg");

		Image image = ImageIO.read(new File("./src/Images/tank3.jpg"));
		//Image image = ImageIO.read(new File("./src/Images/hiresbird.jpg"));

		// Creating buffered image from the given file. NOTE: It's crucial to build the data that way!
		BufferedImage bufferedImage = new BufferedImage(image.getWidth(null), image.getHeight(null), BufferedImage.TYPE_4BYTE_ABGR);
		Graphics g = bufferedImage.getGraphics();
		g.drawImage(image, 0, 0, null);

		System.out.println("BufferedImage: " +bufferedImage.toString());
		
		
		double[][] convert2D = convertTo2DWithoutUsingGetRGB(bufferedImage);
		
		//GaborFilter gfilter = new GaborFilter(5, new double[] {0, Math.PI/2, Math.PI}, 0, 0.5, 1, 3, 3); 
		GaborFilter gfilter = new GaborFilter(5, new double[] {Math.PI}, 0, 0.5, 1, 3, 3); 

		ImageIO.write(gfilter.filter(bufferedImage, null), "jpg", file);
		
		/*
		 * Convolution 2D
		 *  int width, int height, double [][] kernel, int kernelWidth, int kernelHeight
		 */
		
		double[][] result = convolution2D(convert2D, image.getWidth(null), image.getHeight(null), gfilter.getKernelArray(), 3, 3); 
		gfilter.printKernel(result);
		
        BufferedImage filtered = new BufferedImage( image.getWidth(null), image.getHeight(null), BufferedImage.TYPE_4BYTE_ABGR);
        //filtered.setRGB(image.getWidth(null), image.getHeight(null), result);
        int xLength = result.length;
        int yLength =  result[0].length;
        for(int x = 0; x < xLength; x++) {
            for(int y = 0; y < yLength; y++) {
                int rgb = (int)result[x][y]<<16 | (int)result[x][y] << 8 | (int)result[x][y];
                filtered.setRGB(x, y, rgb);
            }
        }
		ImageIO.write(filtered, "jpg", file2);
        System.out.println("\nend");
		//ImageIO.write(new GaborFilter(16, new double[] {0,Math.PI/2, Math.PI}, 0, 0.75, 1, 9, 9).filter(bufferedImage, null), "jpg", file);

	}
}
