function [M, v] = build_M(normals)
    h = size(normals, 1)
    w = size(normals, 2)
    n = h * w
    #M = sparse(2*n, n);
    M = spalloc(2*n, n, n * 4);  # allocate for at least 4*N non-zero values
    v = zeros(2*n, 1);

    for i = 0:(h-1)
        for j = 0:(w-1)
            kr = (i*w + j)*2 + 1;
            kc = i*w + j + 1;
            n = reshape(normals(i+1, j+1, :), 3, 1, 1);
			nx = n(1); ny = n(2); nz = n(3);
            if nz == 0
                nz = 1;
            endif

            # v1
            M(kr, kc) = -nz;
            if j == w - 1
                #M(kr, kc - 1) = nz;
                v(kr, 1) = nx;
            else
                M(kr, kc + 1) = nz;
                v(kr, 1) = -nx;
            endif

            # v2
            M(kr + 1, kc) = -nz;
            if i == h - 1
                #M(kr + 1, kc - w) = nz;
                v(kr + 1, 1) = ny;
            else
                M(kr + 1, kc + w) = nz;
                v(kr + 1, 1) = -ny;
            endif
		endfor
	endfor
end

norms = load("misc/cosas/buda.normals.txt"); norms = reshape(norms, 340, 512, 3);
#h = 2; w = 4; norms = reshape(1:(h*w*3), h, w, 3);
#h = 340; w = 512; norms = rand(h, w, 3);

[M, v] = build_M(norms);
size(M)

A = M'*M;
b = M'*v;
z = A \ b;
#z = zf(1:2:(size(norms,1)*size(norms,2)*2));
z = reshape(z, size(norms, 1), size(norms, 2));
csvwrite("z2.txt", z);
