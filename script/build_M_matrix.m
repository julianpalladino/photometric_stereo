function [M, v] = build_M(normals)
    h = size(normals, 1)
    w = size(normals, 2)
    n = h * w
    M = spalloc(2*n, 2*n, 4*n)  # allocate 4x n for non-zero values
    v = zeros(2*n, 1);

    for i = 0:(h-1)
        for j = 0:(w-1)
            k = (i*w + j)*2 + 1;
            n = reshape(normals(i+1, j+1, :), 3, 1, 1);
			nx = n(1); ny = n(2); nz = n(3);
            if nz == 0
                nz = 10;
            endif

            # v1
            if j == 0
                M(k, k) = -nz;
                M(k, k + 1) = nz;
            else
                M(k, k - 1) = -nz;
                M(k, k + 1) = nz;
			endif
			v(k, 1) = -nx;

            # v2
            if j == 0
                M(k + 1, k) = -nz;
                if j != w - 1
                    M(k + 1, k + 2) = nz;
				endif
            else
                M(k + 1, k - 1) = -nz;
                if j != w - 1
                    M(k + 1, k + 2) = nz;
				endif
			endif
			v(k + 1, 1) = -ny;
		endfor
	endfor
end

norms = load("misc/cosas/buda.normals.txt"); norms = reshape(norms, 340, 512, 3);
#norms = load("pepe.txt"); norms = reshape(norms, 340, 512, 3);

# para recortar
#norms = norms(251:300, 351:400, :);

source("script/filtrar_z.m");

#z = [0 1 5 2 6 3 7 4 8 9 5 6 10 7 11 8 12 9 13 14];
#z = z';

[M, v] = build_M(norms);

A = M'*M;
b = M'*v;
zf = A \ b;
z = filtrar(size(norms, 1), size(norms, 2), zf);
z = reshape(z, size(norms, 1), size(norms, 2));
csvwrite("z.txt", z);
