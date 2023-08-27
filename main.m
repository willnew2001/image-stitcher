clear;

file_1 = 'image1.png';
file_2 = 'image2.png';

% Resize images for speed %
scale_factor = 8;
orig_1 = imresize(double(imread(file_1))/255,(1/scale_factor));
orig_2 = imresize(double(imread(file_2))/255,(1/scale_factor));
gray_1 = rgb2gray(orig_1);
gray_2 = rgb2gray(orig_2);

% RNG %
seed = 1234;
seeded = false;
if seeded
    rng(seed);
end

% Find keypoints on images %
disp("Finding keypoints...");
keypoints_1 = find_keypoints(gray_1);
keypoints_2 = find_keypoints(gray_2);

% Create row and column vectors of all keypoints on images %
[row_1,col_1,~] = find(keypoints_1==1);
[row_2,col_2,~] = find(keypoints_2==1);

% Find best matches for all kp in image 1 %
disp("Finding best matches...");
width = 4; % Region used to describe keypoint (width of 4 gives 9x9 area)
matches_1 = zeros(size(row_1,1),4);
for j=1:size(row_1,1)

    % Collect patch %
    kp = orig_1(row_1(j)-width:row_1(j)+width,col_1(j)-width:col_1(j)+width,:);

    % Compute distances from all keypoints in image 2 %
    dists = zeros(size(row_2,1),1);
    for i=1:size(row_2,1)
        comp = orig_2(row_2(i)-width:row_2(i)+width,col_2(i)-width:col_2(i)+width,:);
        dists(i) = calc_dist(kp,comp);
    end
    
    % Record best match (row,col) pair %
    [sort_dists, I] = sort(dists);
    org = [row_2(I) col_2(I) sort_dists];
    matches_1(j,1) = org(1,1);
    matches_1(j,2) = org(1,2);
    test = find(dists==org(1,3));
    matches_1(j,3) = test(1);
    matches_1(j,4) = org(1,3);
end

% Find best matches for all kp in image 1 %
matches_2 = zeros(size(row_2,1),4);
for j=1:size(row_2,1)

    % Collect patch %
    kp = orig_2(row_2(j)-width:row_2(j)+width,col_2(j)-width:col_2(j)+width,:);

    % Compute distances from all keypoints in image 1 %
    dists = zeros(size(row_1,1),1);
    for i=1:size(row_1,1)
        comp = orig_1(row_1(i)-width:row_1(i)+width,col_1(i)-width:col_1(i)+width,:);
        dists(i) = calc_dist(kp,comp);
    end
    
    % Record best match (row,col) pair %
    [sort_dists, I] = sort(dists);
    org = [row_1(I) col_1(I) sort_dists];
    matches_2(j,1) = org(1,1);
    matches_2(j,2) = org(1,2);
    test = find(dists==org(1,3));
    matches_2(j,3) = test(1);
    matches_2(j,4) = org(1,3);
end

% Find best match pairs %
disp("Finding best pairs...");
dist_thresh = 2; % any match above this distance is filtered out
true_matches = [];
for i=1:size(matches_1,1)
    if (i == matches_2(matches_1(i,3),3))
        if (matches_1(i,4) < dist_thresh)
            true_matches = [true_matches; i matches_1(i,3) matches_1(i,4)];
        end
    end
end

% This value determines what percent of keypoint correspondences are
% necesary to match with the transformation matrix to deem it acceptable
% for stitching %
ratio = (4/7); 

% Find an optimal transformation matrix %
max_transform = [];
max_count = 0;
N = 0;
disp("Finding optimal transformation matrix...");
while (max_count < size(true_matches,1)*ratio && N < 100)
    % Randomly select 4 point correspondences %
    good_matches = randi([1 size(true_matches,1)], 1, 4);
    
    % Create (row,col) matrices for these points on both images %
    good_coords_1 = [];
    good_coords_2 = [];
    for i=1:size(good_matches,2)
        good_coords_1 = [good_coords_1; row_1(true_matches(good_matches(i),1)) col_1(true_matches(good_matches(i),1))];
        good_coords_2 = [good_coords_2; row_2(true_matches(good_matches(i),2)) col_2(true_matches(good_matches(i),2))];
    end

    %%% FINDING THE TRANSFORMATION MATRIX %%%
    % Create A and B matrices %
    A = zeros(size(good_coords_1,1)*2, 9);
    B = zeros(size(good_coords_1,1)*2, 9);
    for rows=1:2:size(A,1)
        % Set x correspondence in A %
        A(rows,1) = good_coords_1((rows+1)/2,2);
        A(rows,2) = good_coords_1((rows+1)/2,1);
        A(rows,3) = 1;
    
        % Set y correspondence in A %
        A(rows+1,4) = good_coords_1((rows+1)/2,2);
        A(rows+1,5) = good_coords_1((rows+1)/2,1);
        A(rows+1,6) = 1;
    
        % Set x correspondence in B %
        B(rows,7) = good_coords_2((rows+1)/2,2) * good_coords_1((rows+1)/2,2);
        B(rows,8) = good_coords_2((rows+1)/2,2) * good_coords_1((rows+1)/2,1);
        B(rows,9) = good_coords_2((rows+1)/2,2);
    
        % Set y correspondence in B %
        B(rows+1,7) = good_coords_2((rows+1)/2,1) * good_coords_1((rows+1)/2,2);
        B(rows+1,8) = good_coords_2((rows+1)/2,1) * good_coords_1((rows+1)/2,1);
        B(rows+1,9) = good_coords_2((rows+1)/2,1);
    end
    
    % Find solution for h %
    [V,D] = eig((transpose(A)*A) + (transpose(B)*B) - (transpose(B)*A) - (transpose(A)*B));
    h = V(:,1);
    %%% TRANSFORMATION MATRIX H HAS BEEN FOUND %%%

    % Count number of point correspondences that match after transformation %
    MATCH_WNDW = round(size(orig_1,2)/22); % Area around point to check for match
    count = 0;
    for i=1:size(true_matches,1)
        % Keypoints and calculated transformation %
        x_1 = col_1(true_matches(i,1));
        y_1 = row_1(true_matches(i,1));
        x_2 = col_2(true_matches(i,2));
        y_2 = row_2(true_matches(i,2));
        x_hat = calc_x_hat(x_1,y_1,h);
        y_hat = calc_y_hat(x_1,y_1,h);

        % If matches ended up within a few pixels of each other, count it %
        if (x_hat > x_2-MATCH_WNDW && x_hat < x_2+MATCH_WNDW) && (y_hat > y_2-MATCH_WNDW && y_hat < y_2+MATCH_WNDW)
            count = count + 1;
        end
    end

    % If new max correspondence match, record transformation matrix %
    if (count > max_count)
        max_count = count;
        max_transform = h;
    end
    N = N + 1;
end

% Set new h to best transformation matrix %
disp("Transformation matrix found.");
h = max_transform;

% X min: top left projected x, bottom left projected x, left edge of base, aka 1 %
min_x = floor(min([calc_x_hat(1, 1, h), calc_x_hat(1, size(orig_1,1), h), 1]));

% X max: top right projected x, bottom right projected x, right edge of base, 756 %
max_x = ceil(max([calc_x_hat(size(orig_1,2), 1, h), calc_x_hat(size(orig_1,2), size(orig_1,1), h), size(orig_1,2)]));

% Y min: top left projected y, top right projected y, top edge of base, aka 1%
min_y = floor(min([calc_y_hat(1,1,h), calc_y_hat(size(orig_1,2),1,h), 1]));

% Y max: bottom left projected y, bottom right projected y, bottom edge of base, 1008 %
max_y = ceil(max([calc_y_hat(1,size(orig_1,1),h), calc_y_hat(size(orig_1,2),size(orig_1,1),h)]));

% Configure canvas size %
canvas = zeros(abs(min_y)+abs(max_y), abs(min_x)+abs(max_x),3);

% Create matrix to transform canvas points back to image 1 %
h_inv = reshape(transpose(inv(transpose(reshape(h,[3,3])))),[9,1]);

% Linear blending feature %
y2 = calc_y_hat(size(orig_1,2),size(orig_1,1),h)+abs(min_y); 
y1 = calc_y_hat(size(orig_1,2),1,h)+abs(min_y);
x2 = calc_x_hat(size(orig_1,2),size(orig_1,1),h)+abs(min_x);
x1 = calc_x_hat(size(orig_1,2),1,h)+abs(min_x);
run = x2 - x1;
rise = y2 - y1;
slope = rise / run;
b = y2 - slope*x2;

% Project the panoramic image %
for rows=1:size(canvas,1) %y
    for cols=1:size(canvas,2) %x

        % Compute the canvas points in both image 1 and image 2 space %
        im1_x = round(calc_x_hat(cols+min_x, rows+min_y, h_inv));
        im1_y = round(calc_y_hat(cols+min_x, rows+min_y, h_inv));
        im2_x = cols+min_x;
        im2_y = rows+min_y;

        % Plot each pixel accordingly %
        on_image_1 = (im1_x > 0 && im1_x <= size(orig_1,2)) && (im1_y > 0 && im1_y <= size(orig_1,1));
        on_image_2 = (im2_x > 0 && im2_x <= size(orig_1,2)) && (im2_y > 0 && im2_y <= size(orig_1,1));

        if on_image_1 && on_image_2
            inner_bound = abs(min_x); % > this
            outer_bound = round((rows - b)/slope); % <= this
            alpha = (cols-inner_bound)/(outer_bound-inner_bound);

            canvas(rows,cols,:) = (1-alpha)*orig_1(im1_y, im1_x,:) + alpha*orig_2(im2_y, im2_x,:);
        elseif (on_image_1)
            canvas(rows,cols,:) = orig_1(im1_y, im1_x,:);
        elseif (on_image_2)
            canvas(rows,cols,:) = orig_2(im2_y, im2_x,:);
        end
       
    end
end
imshow(canvas);


% Create the side by side of each image %
combo = zeros(size(orig_1,1), 2*size(orig_1,2), size(orig_1,3));
for rows = 1:size(combo,1)
    for cols = 1:size(combo,2)
        if (cols <= size(orig_1,2))
           combo(rows,cols,:) = orig_1(rows,cols,:);
       else
            combo(rows,cols,:) = orig_2(rows,cols-size(orig_2,2),:);
        end
    end
end

% Plot point correspondence lines on image %
figure,imshow(combo);
hold on
for i=1:size(good_matches,2)
    p1 = [good_coords_1(i,1),good_coords_1(i,2)];
    p2 = [good_coords_2(i,1),good_coords_2(i,2)+size(orig_2,2)];
    plot([p1(2),p2(2)],[p1(1),p2(1)], 'Color','r','LineWidth',2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x_hat = calc_x_hat(x, y, h)
    x_hat = ((h(1)*x) + (h(2)*y) + h(3)) / ((h(7)*x) + (h(8)*y) + h(9));
end

function y_hat = calc_y_hat(x, y, h)
    y_hat = ((h(4)*x) + (h(5)*y) + h(6)) / ((h(7)*x) + (h(8)*y) + h(9));
end

function dist = calc_dist(A, B)
    dist_sum = 0;
    for i=1:size(A,1)
        for j=1:size(A,2)
            for c=1:size(A,3)
                dist_sum = dist_sum + power(A(i,j,c) - B(i,j,c),2);
            end
        end
    end
    dist = sqrt(dist_sum);
end

function keypoints = find_keypoints(image)
    im1 = image;
    
    % Variables %
    sigma_o = 1.6;
    num_octaves = 4;
    num_scales = 5;
    
    % Define array to hold all scales and octaves %
    smallest_rows = round(size(im1,1)/power(2,num_octaves));
    smallest_cols = round(size(im1,2)/power(2,num_octaves));
    images = zeros(num_octaves,num_scales,smallest_rows*smallest_cols);
    keypoints = zeros(size(im1,1), size(im1,2));
    margin = 4;
    win_size = 1;
    
    % Calculate each scale and octave %
    for octave=1:num_octaves
    
        local_scale = zeros(num_scales, size(im1,1)*size(im1,2));
    
        for scale=1:num_scales
    
            % Calculate sigma %
            sigma = power(2, octave-1) * power(sqrt(2), scale-1) * sigma_o;
    
            % Calculate w %
            if mod(ceil(3*sigma),2) == 1
                w = ceil(3*sigma);
            else
                w = floor(3*sigma);
            end
    
            % Store image in array %
            images(octave,scale,:) = reshape(imresize(imgaussfilt(im1, sigma, FilterSize=w), [smallest_rows, smallest_cols]),[],1);
            local_scale(scale,:) = reshape(imresize(imgaussfilt(im1,sigma,FilterSize=w), [size(im1,1), size(im1,2)]),[],1);
        end
    
        % Compute keypoints %
        for local=1:num_scales
            left = zeros(1,1);
            center = reshape(local_scale(local,:), [size(im1,1), size(im1,2)]);
            right = zeros(1,1);
            if local ~= 1
                left = reshape(local_scale(local-1,:), [size(im1,1), size(im1,2)]);
            end
            if local ~= num_scales
                right = reshape(local_scale(local+1,:), [size(im1,1), size(im1,2)]);
            end
        
            for rows=margin+1:size(center,1)-margin
                for cols=margin+1:size(center,2)-margin
                    if is_local_extrema(center, left, right, rows, cols)
                        keypoints(rows*(power(2,octave-1)),cols*(power(2,octave-1))) = 1;
                    end
                end
            end
        end
        
        % Iteratively subsample %
        im1 = im1(1:2:end,1:2:end);
    end
    
    
    % Filter out edge points %
    edges = edge(image,'canny',0.1);
    for rows=1:size(image,1)
        for cols=1:size(image,2)
            if (keypoints(rows,cols) == 1 && edges(rows,cols) == 1)
                keypoints(rows,cols) = 0;
            end
        end
    end
    
    % Remove keypoints in areas of low contrast %
    for rows=win_size+1:size(image,1)-win_size
        for cols=win_size+1:size(image,2)-win_size
            if keypoints(rows,cols) == 1
                window = image(rows-win_size:rows+win_size, cols-win_size:cols+win_size);
                if (std(reshape(window,[],1)) < 0.08)
                    keypoints(rows,cols) = 0;
                end
            end
        end
    end
end

function is = is_local_extrema(c,l,r, row_coord, col_coord)
    comp = zeros(26,1) + c(row_coord, col_coord);

    % If the center matrix is not scale 1 %
    % Grab left scale pixel values %
    counter = 1;
    if (size(l,1) ~= 1 && size(l,2) ~= 1)
        for rows=row_coord-1:row_coord+1
            for cols=col_coord-1:col_coord+1
                comp(counter) = l(rows,cols);
                counter = counter + 1;
            end
        end
    else
        counter = counter + 9;
    end

    % Grab center scale pixel values %
    for rows=row_coord-1:row_coord+1
        for cols=col_coord-1:col_coord+1
            if ~(rows == row_coord && cols == col_coord)
                comp(counter) = c(rows,cols);
                counter = counter + 1;
            end
        end
    end

    % If the center matrix is not the last scale %
    % Grab right scale pixel values %
    if (size(r,1) ~= 1 && size(r,2) ~= 1)
        for rows=row_coord-1:row_coord+1
            for cols=col_coord-1:col_coord+1
                comp(counter) = r(rows,cols);
                counter = counter + 1;
            end
        end
    end

    max = c(row_coord, col_coord);
    min = c(row_coord, col_coord);

    for i=1:size(comp,1)
        if comp(i) < min
            min = comp(i);
        end
        if comp(i) > max
            max = comp(i);
        end
    end

    if (max == c(row_coord, col_coord) || min == c(row_coord, col_coord))
        is = 1;
    else
        is = 0;
    end
end