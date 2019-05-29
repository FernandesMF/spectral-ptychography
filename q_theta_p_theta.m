% Observáveis q_\theta e p_\theta


theta = 0:30:180;

q_thetas = cell(length(theta),1);

a = annihilfock_mobral(4);

for i=1:length(theta)
            q_thetas{i} = (sqrt(2)/2)*((cosd(theta(i))-i*sind(theta(i)))*a+(cosd(theta(i))+i*sind(theta(i)))*a');
end

            
