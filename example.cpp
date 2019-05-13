#include <stdio.h>
#include <opencv2/opencv.hpp>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <vector>

using namespace cv;
using namespace std;

extern "C"{
#include "vl/sift.h"
};

typedef struct {
    int k1;
    int k2;
    double score;
} Pair;

void detectDescriptor(IplImage *Image, int noctaves, int nlevels, int o_min, vector<float *> &Descriptors, vector<Point2f> &pixels)
{
    // noctaves=(int)(log(min)/log(2));
    vl_sift_pix *ImageData=new vl_sift_pix[Image->height*Image->width];
    unsigned char *Pixel;
    for (int i=0;i<Image->height;i++)
    {
        for (int j=0;j<Image->width;j++)
        {
            Pixel=(unsigned char*)(Image->imageData+i*Image->width+j);
            ImageData[i*Image->width+j]=*(Pixel);
        }
    }
    VlSiftFilt *SiftFilt=NULL;
    SiftFilt=vl_sift_new(Image->width,Image->height,noctaves,nlevels,o_min);

    if (vl_sift_process_first_octave(SiftFilt,ImageData)!=VL_ERR_EOF)
    {
        while (true)
        {
            //计算每组中的关键点
            vl_sift_detect(SiftFilt);
            // iterate all points in the same octave
            VlSiftKeypoint *pKeyPoint=SiftFilt->keys;
            for (int i=0;i<SiftFilt->nkeys;i++)
            {
                VlSiftKeypoint TemptKeyPoint=*pKeyPoint;
                pKeyPoint++;
                //cvDrawCircle(Image,cvPoint(TemptKeyPoint.x,TemptKeyPoint.y),TemptKeyPoint.sigma/2,CV_RGB(255,0,0));
                //idx++;

                double angle[1];
                int angleCount=vl_sift_calc_keypoint_orientations(SiftFilt,angle,&TemptKeyPoint);
                if (angleCount)
                {
                    float *Descriptor=new float[128];
                    vl_sift_calc_keypoint_descriptor(SiftFilt,Descriptor,&TemptKeyPoint,angle[0]);
                    Descriptors.push_back(Descriptor);
                    pixels.push_back(Point2f(TemptKeyPoint.x, TemptKeyPoint.y));
                }
            }
            //下一阶
            if (vl_sift_process_next_octave(SiftFilt)==VL_ERR_EOF)
            {
                break;
            }
        }
    }
    vl_sift_delete(SiftFilt);
}

Pair *match(Pair *pairs, float *descr1, float *descr2, int K1, int K2, int ND, float thresh)
{
    int k1, k2;
    /* Loop over 1st image descr. */
    for (k1 = 0; k1 < K1; ++k1, descr1 += ND ) {
        float best = FLT_MAX;
        float second_best = FLT_MAX;
        int bestk = -1;
        /* Loop over 2nd image descr. and find the 1st and 2nd closest descr. */
        for (k2 = 0; k2 < K2; ++k2, descr2 += ND ) {
            int bin;
            float acc = 0;
            /* Compute the square L2 distance between descriptors */
            for (bin = 0 ; bin < ND ; ++bin) {
                float delta = descr1[bin] - descr2[bin];
                acc += delta*delta;
                if (acc >= second_best)
                    break;
            }
            if (acc < best) {
                second_best = best;
                best = acc;
                bestk = k2;
            }
            else if (acc < second_best) {
                second_best = acc;
            }
        }
        /* Rewind */
        descr2 -= ND*K2;
        /* Record the correspondence if the best descr. passes the ratio test */
        if (thresh * best < second_best && bestk != -1) {
            pairs->k1 = k1;
            pairs->k2 = bestk;
            pairs->score = best;
            pairs++;
        }
    }

    return pairs;
}

int main(int argc, char* argv[])
{
    // load image
    string ImagePath_str1="./lena1.jpg";
    char * ImagePath1=(char*)ImagePath_str1.data();
    IplImage *Image1=cvLoadImage(ImagePath1,0);

    string ImagePath_str2="./lena2.jpg";
    char * ImagePath2=(char*)ImagePath_str2.data();
    IplImage *Image2=cvLoadImage(ImagePath2,0);

    // setup timer
    struct  timeval start;
    struct  timeval end;
    unsigned long diff;
    gettimeofday(&start,NULL);

    // set descriptor dimension
    int ND = 128;
    vector<float *> Descriptors1, Descriptors2;
    vector<Point2f> pixels1, pixels2;

    // extract (paremeters can be changed manually)
    detectDescriptor(Image1, 4, 2, 0, Descriptors1, pixels1);
    detectDescriptor(Image2, 4, 2, 0, Descriptors2, pixels2);
    // convert vector<float *> to a raw float list
    int K1 = Descriptors1.size();
    int K2 = Descriptors2.size();
    float *desc1 = new float[K1 * ND];
    float *desc2 = new float[K2 * ND];
    for (int i = 0; i < K1; ++i)
    {
        for (int j = 0; j < ND; ++j)
        {
            *desc1++ = Descriptors1[i][j];
        }

    }
    // Don't forgect to move back!!!
    desc1 -= K1 * ND;
    for (int i = 0; i < K2; ++i)
    {
        for (int j = 0; j < ND; ++j)
        {
            *desc2++ = Descriptors2[i][j];
        }

    }
    // Don't forgect to move back!!!
    desc2 -= K2 * ND;

    // perform match
    Pair* pairs_begin = (Pair*) malloc(sizeof(Pair) * (K1+K2));
    Pair* pairs_iterator = pairs_begin;
    pairs_iterator = match(pairs_iterator, desc1, desc2, K1, K2, ND, 1.5);
    Pair* pairs_end = pairs_iterator;
    // Init match list and score list
    int *out1 = new int[pairs_end - pairs_begin];
    int *out2 = new int[pairs_end - pairs_begin];
    int *scores = new int[pairs_end - pairs_begin];

    // for visualize
    cv::Mat image1, image2;
    image1 = cv::imread("./lena1.jpg");
    image2 = cv::imread("./lena2.jpg");
    cv::hconcat(image1, image2, image2);
    for (pairs_iterator = pairs_begin; pairs_iterator < pairs_end; ++pairs_iterator)
    {
      // matched index: k1 and k2
      *out1++ = pairs_iterator->k1;
      *out2++ = pairs_iterator->k2;
      *scores++ = pairs_iterator->score;
      //cout<<"match: "<<pairs_iterator->k1<<"-"<<pairs_iterator->k2<<endl;
      //cout<<"coord: "<<pixels1[pairs_iterator->k1]<<"-"<<pixels2[pairs_iterator->k2]<<endl;

      // for visualize
      Point2f move256pixel(pixels2[pairs_iterator->k2].x+256, pixels2[pairs_iterator->k2].y);
      cv::line(image2, pixels1[pairs_iterator->k1], move256pixel, cv::Scalar(255, 0, 0), 1.5);
      cv::circle(image2, pixels1[pairs_iterator->k1], 2, cv::Scalar(0, 0, 255), 2);
      cv::circle(image2, move256pixel, 2, cv::Scalar(0, 255, 0), 2);

    }

    // Timer
    gettimeofday(&end,NULL);
    diff = 1000000 * (end.tv_sec-start.tv_sec)+ end.tv_usec-start.tv_usec;
    printf("Time cost: %ld us\n",diff);

    // for visualize
    cv::imshow("demo", image2);
    cv::waitKey(99999);
    return 0;
}
