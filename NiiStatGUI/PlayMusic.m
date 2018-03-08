function PlayMusic
    [y, Fs] = audioread('backgroundMusic_01.mp3'); 
    player = audioplayer(y, Fs);
    play(player)
end