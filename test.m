clc; close all; clear ;

EbN0_vec          =  -10 : 10 : 30;

N_Packets = 25;

N_data = 500;

Bit_length =16;

Symbol_length= Bit_length/2;

for Active_Antenna = 1 : 1 :  8    
    % Set Active Antennas
    
Amplitude_tx1 = ( Active_Antenna > 0 ) * 1;
Amplitude_tx2 = ( Active_Antenna > 1 ) * 1;
Amplitude_tx3 = ( Active_Antenna > 2 ) * 1;
Amplitude_tx4 = ( Active_Antenna > 3 ) * 1;
Amplitude_tx5 = ( Active_Antenna > 4 ) * 1;
Amplitude_tx6 = ( Active_Antenna > 5 ) * 1;
Amplitude_tx7 = ( Active_Antenna > 6 ) * 1;
Amplitude_tx8 = ( Active_Antenna > 7 ) * 1;

BER_EbN0 = ones(length(EbN0_vec),1);               % Memory Allocation

for j = 1 : 1 : length(EbN0_vec)                    % LOOP: EbN0
    EbN0_lin = 10^(EbN0_vec(j)/10);             % Calculate Noise Power
    N0 = 1 / EbN0_lin;
    sigma_n = sqrt(N0);

    BER_packets = ones(length(N_Packets),1);
    for packets = 1:1:N_Packets

        Tbits = randi([0,1],[Bit_length*N_data,1]);
        Transmit_symbols = qammod(Tbits,4,'InputType','bit','UnitAveragePower',true);

        for bits=0:1:N_data-1
            %% ENCODE

            Enc = zeros(Symbol_length*2-1,Symbol_length);

            Tsymbols = Transmit_symbols((bits*Symbol_length+1):((bits+1)*Symbol_length),:);

            A=fft(eye(Symbol_length));

            for Symbol=1:1:Symbol_length
                Enc(:,Symbol)=[zeros(Symbol-1,1);Tsymbols;zeros(Symbol_length-Symbol,1)];
            end

            X = Enc*A;

            c_d = zeros(Symbol_length,1);
            
            X1=zeros(Symbol_length*2-1,Symbol_length);
            
            for t=1:1:Symbol_length
            c_d(t) = sqrt(mean(abs(X(:,t)).^2));
            X1(:,t) =  X(:,t) / c_d(t);                 
            end        

            %% CHANNEL 
            H_channel = sqrt(1/2)*(randn(Symbol_length,1)+1j*(randn(Symbol_length,1)));
            
            H_channel(1) = H_channel(1) * Amplitude_tx1;
            H_channel(2) = H_channel(2) * Amplitude_tx2;
            H_channel(3) = H_channel(3) * Amplitude_tx3;
            H_channel(4) = H_channel(4) * Amplitude_tx4;
            H_channel(5) = H_channel(5) * Amplitude_tx5;
            H_channel(6) = H_channel(6) * Amplitude_tx6;
            H_channel(7) = H_channel(7) * Amplitude_tx7;
            H_channel(8) = H_channel(8) * Amplitude_tx8;
            
            y1((bits*(Symbol_length*2-1)+1):((bits+1)*(Symbol_length*2-1)),:) = X1*H_channel;


            %% NOISE 
            noise = sigma_n/sqrt(2)*(randn(size(y1))+1j*randn(size(y1)));
            Ps = mean(abs(y1).^2);
            Pn = mean(abs(noise).^2);
            Transmit_Signal = y1 + noise;

            %% DECODE
            H_normalised = H_channel./c_d;

            H_bar = A * H_normalised;
            H_Eq = zeros(Symbol_length*2-1,Symbol_length);
            for Symbol=1:1:Symbol_length
            H_Eq(:,Symbol)=[zeros(Symbol-1,1);H_bar;zeros(Symbol_length-Symbol,1)];
            end
            y = Transmit_Signal((bits*15+1):((bits+1)*15),:);
            Rsymbols((bits*8+1):((bits+1)*8),:) = ((H_Eq.')*H_Eq)\H_Eq.'*y;

            test1 = H_Eq* Tsymbols  ;
            test2 = X1*H_channel;
%             plot(test1,test2);                               %DEBUG
         end
         Rbits = qamdemod(Rsymbols,4,'OutputType','bit','UnitAveragePower',true);
         BER_packets(packets) = sum(abs(Tbits- Rbits))/length(Tbits);
    end
    BER_EbN0(j) = mean(BER_packets);
end
figure(2)
semilogy(EbN0_vec,BER_EbN0,'-x','LineWidth',2,'MarkerSize',4);
hold on;
grid()
end