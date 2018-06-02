function [M, v] = build_M(normals, mask)
    [h, w] = size(mask);

    mask_cosa = and(normals(:,:,3) != 0, mask != 0);
    mask_idx = find(mask_cosa');
    [mask_x, mask_y] = find(mask_cosa');
    mask_xy = [mask_x mask_y];

    n = size(mask_idx, 1);

    mask_dicc = uint64(mask');
    mask_dicc(mask_idx) = 1:size(mask_idx, 1);

    #h = size(normals, 1)
    #w = size(normals, 2)
    #n = h * w
    M = sparse(2*n, n);  # allocate for at least 4*N non-zero values
    v = zeros(2*n, 1);

    for i = 1:n
        x = mask_x(i);
        y = mask_y(i);
        #disp([x y]);

        normal = reshape(normals(y, x, :), 3, 1, 1);
        nn = norm(reshape(normals(y, x, :), 3, 1, 1));
        nx = normals(y, x, 1) / nn;
        ny = normals(y, x, 2) / nn;
        nz = normals(y, x, 3) / nn;

        # ecuacion 1
        M(i, i) = -nz;
        if x < w && mask(y, x+1) > 0
            M(i, i+1) = nz;
        endif
        v(i, 1) = -nx;

        # ecuacion 2
        j = i + n;
        M(j, i) = -nz;
        if y < h && mask(y+1, x) > 0
            col = mask_dicc(x, y+1);
            M(j, col) = nz;
        endif
        v(j, 1) = -ny;
    endfor
end

function [z] = build_z(zf, normals, mask)
    mask_cosa = and(normals(:,:,3), mask);
    mask_idx = find(mask_cosa');
    [mask_x, mask_y] = find(mask_cosa');
    mask_xy = [mask_x mask_y];

    n = size(mask_idx, 1);

    [h, w] = size(mask);
    z = zeros(h, w);

    for i = 1:n
        x = mask_x(i);
        y = mask_y(i);
        z(y, x) = zf(i);
    endfor
end

#norms = load("misc/cosas/buda.normals.txt"); norms = reshape(norms, 340, 512, 3);
nx = csvread("buda.nx.txt"); nx = nx(:, 1:end-1);
ny = csvread("buda.ny.txt"); ny = ny(:, 1:end-1);
nz = csvread("buda.nz.txt"); nz = nz(:, 1:end-1);
norms = cat(3, nx, ny, nz);

#mask = imread("misc/buda/buda.mask.ppm");
#[h, w, ~] = size(mask);
#[htrim, wtrim, ~] = size(nx);
#mask = mask((h/2 - htrim/2 + 1):(h/2 + htrim/2), (w/2 - wtrim/2 + 1):(w/2 + wtrim/2));

[M, v] = build_M(norms, mask);
size(M)

#Mposta = M; vposta = v;
#load("M.dat"); M = sparse(M);
#load("v.dat");

A = M'*M;
b = M'*v;

#Aposta = A; bposta = b;
#load("A.dat"); A = sparse(A);
#load("b.dat");

zf = A \ b;
z = build_z(zf, norms, mask);
z(z==0) = min(min(z));
csvwrite("z.txt", z);
