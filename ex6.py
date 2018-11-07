##################################################################
# FILE : ex6.py
# WRITER : Lior Paz, lioraryepaz, 206240996
# EXERCISE : intro2cs ex6 2017-2018
# DESCRIPTION : Cutting edge Image Processing script - you give at an image,
#  it gives you back an adjusted one in the proper angle :) also, if you
# wish, you can use it to apply different filters on the image, and reducing
#  it size!!!! how cool!!!
##################################################################

import ex6_helper as ph
# let us use helping functions from cs team
import sys
# let us use cmd arguments
from collections import Counter
# for the counting in Otsu
import copy
# for list copying
import math

EDGES_FILTER = [[-1 / 8, -1 / 8, -1 / 8], [-1 / 8, 1, -1 / 8],
                [-1 / 8, -1 / 8, -1 / 8]]

NUMBER_OF_ARGUMENTS = 3
# for sys_arg

ERROR_MESSAGE = 'ERROR: invalid number of parameters. Please enter: ' \
                'ex6.py<image_source><output><max_diagonal>'

BLACK = 0

WHITE = 255

RIGHT_ANGLE = 90


def otsu(image):
    """Gives us the optimal threshold for the image"""
    count_dict = Counter(x for sublist in image for x in sublist)
    # creates a dictionary which contains the counting for every
    # pixel type
    num_black_pix = 0
    # in black we start from 0
    sum_black_pix = 0
    num_white_pix = sum(count_dict.values())
    # in white we begin with all of the values
    sum_white_pix = 0
    for val in count_dict:
        sum_white_pix += val * count_dict[val]
    intra_variance_dict = {0: 0}
    for threshold in range(1, 256):
        # for every threshold, im simply 'moving' the current pixel type from
        # the 'white side' to the 'black side'
        num_black_pix += count_dict[threshold - 1]
        sum_black_pix += (threshold - 1) * count_dict[threshold - 1]
        num_white_pix -= count_dict[threshold - 1]
        sum_white_pix -= (threshold - 1) * count_dict[threshold - 1]
        if num_black_pix == 0:
            mean_black = 0
        else:
            mean_black = sum_black_pix / num_black_pix
        if num_white_pix == 0:
            mean_white = 0
        else:
            mean_white = sum_white_pix / num_white_pix
        intra_variance = num_black_pix * num_white_pix * (
                mean_black - mean_white) ** 2
        intra_variance_dict.update({threshold: intra_variance})
        # the dict that contains the intra val for every pixel type
    # the max val is the one we want
    return max(intra_variance_dict, key=intra_variance_dict.get)


def threshold_filter(image):
    """Gives us the image in solid black & white colors according to otsu
    threshold"""
    threshold_image = copy.deepcopy(image)
    # creating a whole new image list so we wont change the original
    ot = otsu(image)
    # gives us the otsu val for the image
    for row in range(len(threshold_image)):
        for column in range(len(threshold_image[0])):
            if threshold_image[row][column] < ot:
                threshold_image[row][column] = BLACK
            else:
                threshold_image[row][column] = WHITE
    return threshold_image


def calc(filter, u_l, u_c, u_r, c_l, c, c_r, d_l, d_c, d_r):
    """Performing the calculations for apply_filter"""
    temp = u_l * filter[0][0] + u_c * filter[0][1] \
           + u_r * filter[0][2] + c_l * filter[1][0] + \
           c * filter[1][1] + c_r * filter[1][2] + \
           d_l * filter[2][0] + d_c * filter[2][1] + \
           d_r * filter[2][2]
    # alter the value correctly
    return min(int(math.fabs(temp)), WHITE)


def detect_edges(image):
    """Is using an algorithm to find edges"""
    edges_image = apply_filter(image, EDGES_FILTER)
    # the edges filter is equal actually to the requirement - (pixel - average
    # of 8 neighbours)
    return edges_image


def apply_filter(image, filter):
    """Alter every pixel in the image in a specific way, according to its
    neighbours"""
    # important - this func code is very long for a reason - i wanted to
    # prevent my 'big' for loop from going over un-necessary ifs, because
    # there is most of my running-time, so i created different use-cases
    new_image = copy.deepcopy(image)
    # creating a whole new image list so we wont change the original
    end_row = len(image) - 1
    end_column = len(image[0]) - 1
    for row in range(1, end_row):
        # going over the first column without edges
        c = image[row][0]
        u_l = c
        u_c = image[row - 1][0]
        u_r = image[row - 1][1]
        c_l = c
        c_r = image[row][1]
        d_l = c
        d_c = image[row + 1][0]
        d_r = image[row + 1][1]
        new_image[row][0] = calc(filter, u_l, u_c, u_r, c_l, c, c_r, d_l, d_c,
                                 d_r)
        # going over the last column without edges
        c = image[row][end_column]
        u_l = image[row - 1][end_column - 1]
        u_c = image[row - 1][end_column]
        u_r = c
        c_l = image[row][end_column - 1]
        c_r = c
        d_l = image[row + 1][end_column - 1]
        d_c = image[row + 1][end_column]
        d_r = c
        new_image[row][end_column] = calc(filter, u_l, u_c, u_r, c_l, c,
                                          c_r, d_l, d_c, d_r)
        # this is the Big loop - going over the rest of the columns without
        # edges
        for column in range(1, end_column):
            u_l = image[row - 1][column - 1]
            u_c = image[row - 1][column]
            u_r = image[row - 1][column + 1]
            c_l = image[row][column - 1]
            c = image[row][column]
            c_r = image[row][column + 1]
            d_l = image[row + 1][column - 1]
            d_c = image[row + 1][column]
            d_r = image[row + 1][column + 1]
            new_image[row][column] = calc(filter, u_l, u_c, u_r, c_l,
                                          c, c_r, d_l, d_c, d_r)
    for column in range(1, end_column):
        # going over the first row without edges
        c = image[0][column]
        u_l = c
        u_c = c
        u_r = c
        c_l = image[0][column - 1]
        c_r = image[0][column + 1]
        d_l = image[1][column - 1]
        d_c = image[1][column]
        d_r = image[1][column + 1]
        new_image[0][column] = calc(filter, u_l, u_c, u_r, c_l, c, c_r, d_l,
                                    d_c, d_r)
        # going over the last row without edges
        c = image[end_row][column]
        u_l = image[end_row - 1][column - 1]
        u_c = image[end_row - 1][column]
        u_r = image[end_row - 1][column + 1]
        c_l = image[end_row][column - 1]
        c_r = image[end_row][column + 1]
        d_l = c
        d_c = c
        d_r = c
        new_image[end_row][column] = calc(filter, u_l, u_c, u_r, c_l,
                                          c, c_r, d_l, d_c, d_r)
    # changing the top-left pixel
    c = image[0][0]
    u_l = c
    u_c = c
    u_r = c
    c_l = c
    c_r = image[0][1]
    d_l = c
    d_c = image[1][0]
    d_r = image[1][1]
    new_image[0][0] = calc(filter, u_l, u_c, u_r, c_l, c, c_r, d_l, d_c, d_r)
    # changing the bottom-right pixel
    c = image[end_row][end_column]
    u_l = image[end_row - 1][end_column - 1]
    u_c = image[end_row - 1][end_column]
    u_r = c
    c_l = image[end_row][end_column - 1]
    c_r = c
    d_l = c
    d_c = c
    d_r = c
    new_image[end_row][end_column] = calc(filter, u_l, u_c, u_r, c_l, c,
                                          c_r, d_l, d_c, d_r)
    # changing the top-right pixel
    c = image[0][end_column]
    u_l = c
    u_c = c
    u_r = c
    c_l = image[0][end_column - 1]
    c_r = c
    d_l = image[1][end_column - 1]
    d_c = image[1][end_column]
    d_r = c
    new_image[0][end_column] = calc(filter, u_l, u_c, u_r, c_l, c,
                                    c_r, d_l, d_c, d_r)
    # changing the bottom-left pixel
    c = image[end_row][0]
    u_l = c
    u_c = image[end_row - 1][0]
    u_r = image[end_row - 1][1]
    c_l = c
    c_r = image[end_row][1]
    d_l = c
    d_c = c
    d_r = c
    new_image[end_row][0] = calc(filter, u_l, u_c, u_r, c_l, c,
                                 c_r, d_l, d_c, d_r)
    return new_image


def downsample_by_3(image):
    """Reducing the resolution of an image times 3"""
    down_tri_image = []
    # jumps by 3 to the center of the next new pixel, and ending in '-1' to
    # prevent crashes if the image isn't 3^n*3^n
    for row in range(1, len(image) - 1, 3):
        temp = []
        for column in range(1, len(image[0]) - 1, 3):
            # gives me the average of the 9 pixels - the new value of the
            # new pixel
            temp.append(int(((image[row][column] + image[row - 1][column - 1]
                              + image[row - 1][column] + image[row - 1][
                                  column + 1] +
                              image[row][column - 1] + image[row][column + 1] +
                              image[row + 1][column - 1] + image[row + 1][
                                  column] +
                              image[row + 1][column + 1]) / 9)))
        down_tri_image.append(temp)
    return down_tri_image


def downsample(image, max_diagonal_size):
    """Reducing the resolution of an image to a requested size"""
    # i used here recursion - im  using downsample_by_3 func, as base case
    # of when the image is in the requested size
    if math.sqrt(math.pow(len(image), 2) + math.pow(len(image[0]),2)) <= \
            max_diagonal_size:
        return image
    else:
        return downsample(downsample_by_3(image), max_diagonal_size)


def current_rank(image, angle, distance, top):
    """Calculate the rank for a specific line"""
    line_index = ph.pixels_on_line(image, angle * math.pi / 180, distance, top)
    # calling the help func in order to get the specific pixels i need
    white_lst = []
    temp_dist = 0
    current_rank = 0
    for pixel in line_index:
        # creating a renwed list with only white pixels in it
        if image[pixel[0]][pixel[1]] == WHITE:
            white_lst.append(pixel)
    for p in range(len(white_lst) - 1):
        # here we calculate the distance between every following white
        # pixels - as long as they are close enough to each other, we keep
        # adding the small distances to one large distance, of every segment
        #  in separate. when we reach a 'space', we begin new segment
        d = math.sqrt(math.pow((white_lst[p][0] - white_lst[p + 1][0]),
                               2) + math.pow(
            white_lst[p][1] - white_lst[p + 1][1], 2))
        if d <= 2:
            temp_dist += d
        else:
            current_rank += math.pow(temp_dist, 2)
            temp_dist = 0
    # this command line is for the end of every line the function is working on
    current_rank += math.pow(temp_dist, 2)
    return current_rank


def get_angle(image):
    """The func gives us the dominant angle in a specific image, according
    to a special formula we developed"""
    diagonal = math.sqrt(math.pow(len(image), 2) + math.pow(len(image[0]), 2))
    # Pythagoras
    angle_dict = {}
    # now we begin to use different use cases - because for 1-89 we have 2
    # lines, and for distance 0 we have special requirements
    # angle == 0
    angle_rank = 0
    for distance in range(int(diagonal) + 1):
        angle_rank += current_rank(image, 0, distance, True)
    angle_dict.update({0: angle_rank})
    # angle is 1-89
    for angle in range(1, RIGHT_ANGLE):
        angle_rank = 0
        # distance == 0
        angle_rank += current_rank(image, angle, 0, True)
        # distance > 0
        for distance in range(1, int(diagonal) + 1):
            angle_rank += current_rank(image, angle, distance, True)
            angle_rank += current_rank(image, angle, distance, False)
        angle_dict.update({angle: angle_rank})
    # angle == 90
    angle_rank = 0
    for distance in range(int(diagonal) + 1):
        angle_rank += current_rank(image, RIGHT_ANGLE, distance, True)
    angle_dict.update({RIGHT_ANGLE: angle_rank})
    # angle == 91-180
    for angle in range(91, 180):
        angle_rank = 0
        for distance in range(1, int(diagonal) + 1):
            angle_rank += current_rank(image, angle, distance, True)
        angle_dict.update({angle: angle_rank})
    # well give the angle with the highest rank
    return max(angle_dict, key=angle_dict.get)


def rotate(image, angle):
    """The func receive an image and rotates it in a desired angle"""
    radian = math.radians(angle)
    radian_b = math.radians(angle - RIGHT_ANGLE)
    sin_angle = math.sin(radian)
    cos_angle = math.cos(radian)
    length = len(image)
    length_0 = len(image[0])
    # gives us the different rotated frames in bigger size
    if angle < RIGHT_ANGLE:
        side_a = length_0 * cos_angle + length * sin_angle
        side_b = length_0 * sin_angle + length * cos_angle
    if angle == RIGHT_ANGLE:
        side_a = length
        side_b = length_0
    if angle > RIGHT_ANGLE:
        side_a = length * math.cos(radian_b) + length_0 * math.sin(radian_b)
        side_b = length * math.sin(radian_b) + length_0 * math.cos(radian_b)
    rotated_image = []
    # gives me a start point of black pixels
    for row in range(int(math.fabs(side_b))):
        temp = []
        for pixel in range(int(math.fabs(side_a))):
            temp.append(BLACK)
        rotated_image.append(temp)
    x = len(rotated_image[0]) / 2
    y = len(rotated_image) / 2
    old_x = length_0 / 2
    old_y = length / 2
    for row in range(len(rotated_image)):
        for pixel in range(len(rotated_image[0])):
            # for every pixel - if its in the range of the original image
            # and not supposed to be black 'edge' - we use simple
            # trigonometry to find the origin index
            new_row = row - y
            new_pixel = pixel - x
            old_row = new_pixel * sin_angle + new_row * cos_angle
            old_pixel = new_pixel * cos_angle - new_row * sin_angle
            old_row_centered = int(old_row + old_y)
            old_pixel_centered = int(old_pixel + old_x)
            if (-1 < old_row_centered < old_y * 2) and (-1 <
                                            old_pixel_centered < old_x * 2):
                rotated_image[row][pixel] = image[old_row_centered][
                    old_pixel_centered]
    return rotated_image


def make_correction(image, max_diagonal):
    """Perform a complete correction to an image"""
    angle = get_angle(threshold_filter(detect_edges(threshold_filter(
        downsample(image, max_diagonal)))))
    return rotate(image, angle)


if __name__ == "__main__":
    if len(sys.argv) == NUMBER_OF_ARGUMENTS + 1:
        # the if check for wrong number of arguments
        image_source = sys.argv[1]
        output = sys.argv[2]
        max_diagonal = int(sys.argv[3])
        image = ph.load_image(image_source)
        new_image = make_correction(image, max_diagonal)
        ph.save(new_image, output)
    else:
        print(ERROR_MESSAGE)
