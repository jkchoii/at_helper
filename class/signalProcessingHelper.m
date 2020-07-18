classdef signalProcessingHelper
    % Acoustics Toolbox�� ���ǵ� ��ȣ�� �ٷ� �� ������ �Ǵ� Ŭ����
    % methods:
    %   public:
    %       getCovarianceMatrix(*) - ���л� ��� ����
    %       getAverageNoisySnapshot(*) - SNR ����þ� ������ ���� ��ȣ L�� ��ճ�
    %       getGaussianNoise(*) - SNR ����þ� ���� ��
    
    methods
        
        function [STFT, f, t] = stft(obj, x, win, hop, nfft, fs, varargin)

        % function: [STFT, f, t] = stft(x, win, hop, nfft, fs)
        %
        % Input:
        % x - signal in the time domain
        % win - analysis window function
        % hop - hop size
        % nfft - number of FFT points
        % fs - sampling frequency, Hz
        %
        % Output:
        % STFT - STFT-matrix (only unique points, time 
        %        across columns, frequency across rows)
        % f - frequency vector, Hz
        % t - time vector, s

        % representation of the signal as column-vector
        x = x(:);

        % determination of the signal length 
        xlen = length(x);

        % determination of the window length
        wlen = length(win);

        % stft matrix size estimation and preallocation
        if ~isempty(varargin)
            if strcmpi(varargin{1}, 'SYMMETRIC') 
                NUP = nfft;
            else
                error('�Ű����� Ȯ��')
            end
        else
            NUP = ceil((1+nfft)/2);     % calculate the number of unique fft points
            
        end
        
        L = 1+fix((xlen-wlen)/hop); % calculate the number of signal frames
        STFT = zeros(NUP, L);       % preallocate the stft matrix

        % STFT (via time-localized FFT)
        for l = 0:L-1
            % windowing
            xw = x(1+l*hop : wlen+l*hop).*win;

            % FFT
            X = fft(xw, nfft);

            % update of the stft matrix
            STFT(:, 1+l) = X(1:NUP);
        end

        % calculation of the time and frequency vectors
        t = (wlen/2:hop:wlen/2+(L-1)*hop)/fs;
        f = (0:NUP-1)*fs/nfft;

        end
        
        function cov = getCovarianceMatrix(obj, p)
            % ���л� ��� �����
            % (����) ���л� ����� ���� ������ ������ �߰��ϸ� ��Ī ����� �ƴ�
            % ������ �߰��� ��ȣ ���͸� �Է����� �־�� ��
            % p: ������ �߰��� ��ȣ ����
            % cov: ���� ��Ī ���
            
            assert(isvector(p), '���Ͱ� �ƴ�');
            
            p = p ./ max(p);
            cov = p * p';

        end
        
        function new_s = getAverageNoisySnapshots(obj, s, L, SNR)
            %METHOD1 �� �޼����� ��� ���� ��ġ
            %   �ڼ��� ���� ��ġ
            
            assert(isvector(s), '���Ͱ� �ƴ�');
            new_s = zeros(length(s), 1);
            for idx = 1:L
                new_s = new_s + getGaussianNoise(obj, s, SNR);
            end
            new_s = new_s / L;
            new_s = new_s ./ max(new_s);

        end
        
        function new_s = getGaussianNoise(obj, s, SNR)
            %METHOD1 �� �޼����� ��� ���� ��ġ
            %   �ڼ��� ���� ��ġ
            
            assert(isvector(s), '���Ͱ� �ƴ�');
            new_s = awgn(s, SNR, 'measured');
            new_s = new_s ./ max(new_s);
        end
        
    end
end
