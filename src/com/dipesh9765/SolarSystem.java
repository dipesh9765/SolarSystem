package com.dipesh9765;

import javax.swing.*;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.geom.*;
import java.util.Vector;


interface Host {

    double getAbsPos(byte mAxis);

    double getMass();
}


class Coord {

    public static final byte X = 1, Y = 0;
}

public class SolarSystem {

    private static final double dt = 0.001;


    public SolarSystem() {
    }


    public static void main(String[] arguments) {
        SkyFrame skyFrame;
        int numOfPlanets = 4;

        if (arguments.length == 0) skyFrame = new SkyFrame(numOfPlanets, dt);
        else {
            if (arguments[0].equals("h")) {
                usage();
                System.exit(1);
            }
            try {
                numOfPlanets = Integer.parseInt(arguments[0]);
                numOfPlanets = Math.abs(numOfPlanets);
                skyFrame = new SkyFrame(numOfPlanets, 0);
            } catch (Exception e) {
                System.err.println(e);
                usage();
                System.exit(0);
                skyFrame = new SkyFrame(0, 0);
            }
        }

        while (true) {
            skyFrame.getPanel().repaint();
        }
    }


    public static Color getColor(int i) {
        Color color;
        int j = i % 8;

        color = switch (j) {
            case 0 -> Color.red;
            case 1 -> Color.blue;
            case 2 -> Color.green;
            case 3 -> Color.pink;
            case 4 -> Color.orange;
            case 5 -> Color.cyan;
            case 6 -> Color.magenta;
            case 7 -> Color.lightGray;
            default -> Color.gray;
        };
        return color;
    }


    private static void usage() {
        String explaination = """
                This program displays a simulated solar system in 2 dimensions.
                Buttons along the top of the GUI allow the user to show or hide
                various physical properties of the planets. These properties are
                the acceleration and velocity vectors, the percent kinetic and potential
                energy, and the trajectory the planet has taken. In addition, two
                buttons allow the user to zoom in or out for better viewing. Finally,
                there is a sliding scale on the right of the GUI that allows the user
                to alter the mass of the sun.""";

        System.out.println("Usage :\n" + "java SolarSystem [options] <number of planets>\n" + "\n" + "options :\n" + "h\t\tprint this message\n\n" + explaination);
    }
}


class SkyFrame extends JFrame implements ChangeListener, ActionListener {

    private static final double maxSunMassMultiplier = 1.5;
    private static final double minSunMassMultiplier = 0.5;

    private static final int xPixels = 1000;
    private static final int yPixels = 1000;

    private static final int maxScaleFactor = 5;

    private static int lengthScaleFactor = 1;
    private final String zoomInString = "Zoom In", zoomOutString = "Zoom Out";
    private final String showAcc = "Show Acceleration";
    private final String showVel = "Show Velocity";
    private final String showEnergy = "Show Energy";
    private final String showPos = "Trace Position";
    private final String freeze = "Freeze";
    private final SkyPanel skyPanel;
    private final Sun sun;
    private final JSlider slider;
    private final JButton zoomIn;
    private final JButton zoomOut;
    private final JButton accButton;
    private final JButton velButton;
    private final JButton energyButton;
    private final JButton positionButton;
    private final JButton freezeButton;

    public SkyFrame(int mNumberOfPlanets, double dt) {
        super("Solar System");

        sun = new Sun(75, 1000, xPixels / 2, yPixels / 2);

        setSize(xPixels, yPixels);
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

        skyPanel = new SkyPanel(dt, xPixels, mNumberOfPlanets, sun);
        skyPanel.setBackground(Color.black);

        BorderLayout MasterPanelLayout = new BorderLayout();
        Container pane = getContentPane();
        pane.setLayout(MasterPanelLayout);
        pane.add(skyPanel, BorderLayout.CENTER);

        slider = new JSlider(JSlider.VERTICAL, (int) (sun.getMass() * minSunMassMultiplier), (int) (sun.getMass() * maxSunMassMultiplier), (int) sun.getMass());
        slider.setMajorTickSpacing((int) (sun.getMass() / 20));
        slider.setMinorTickSpacing((int) (sun.getMass() / 40));
        slider.setPaintTicks(true);
        slider.setPaintLabels(true);
        slider.addChangeListener(this);
        pane.add(slider, BorderLayout.EAST);

        JPanel buttonPanel = new JPanel();
        zoomIn = new JButton(zoomInString);
        zoomIn.addActionListener(this);
        zoomIn.setEnabled(false);
        buttonPanel.add(zoomIn);

        zoomOut = new JButton(zoomOutString);
        zoomOut.addActionListener(this);
        zoomOut.setEnabled(true);
        buttonPanel.add(zoomOut);

        accButton = new JButton(showAcc);
        accButton.addActionListener(this);
        buttonPanel.add(accButton);

        velButton = new JButton(showVel);
        velButton.addActionListener(this);
        buttonPanel.add(velButton);

        energyButton = new JButton(showEnergy);
        energyButton.addActionListener(this);
        buttonPanel.add(energyButton);

        positionButton = new JButton(showPos);
        positionButton.addActionListener(this);
        buttonPanel.add(positionButton);

        freezeButton = new JButton(freeze);
        freezeButton.addActionListener(this);
        buttonPanel.add(freezeButton);

        pane.add(buttonPanel, BorderLayout.NORTH);

        setContentPane(pane);
        setVisible(true);
    }


    public static int getScaleFactor() {
        return lengthScaleFactor;
    }


    public void actionPerformed(ActionEvent event) {
        String command = event.getActionCommand();

        String hideAcc = "Hide Acceleration";
        String hideVel = "Hide Velocity";
        String hideEnergy = "Hide Energy";
        String hidePos = "Hide Trace";
        String unFreeze = "Unfreeze";
        if (command.equals(zoomOutString)) {
            zoomIn.setEnabled(true);
            ++lengthScaleFactor;
            energyButton.setEnabled(false);
            skyPanel.displayEnergy(false);
            skyPanel.populateStarGroup(lengthScaleFactor);
            if (lengthScaleFactor > maxScaleFactor) zoomOut.setEnabled(false);
        } else if (command.equals(zoomInString)) {
            --lengthScaleFactor;
            if (lengthScaleFactor == 1) {
                zoomIn.setEnabled(false);
                energyButton.setEnabled(true);
                if (energyButton.getText().equals(hideEnergy)) skyPanel.displayEnergy(true);
            }
            zoomOut.setEnabled(true);
        } else if (command.equals(showAcc)) {
            skyPanel.displayAcc(true);
            accButton.setText(hideAcc);
        } else if (command.equals(showVel)) {
            skyPanel.displayVel(true);
            velButton.setText(hideVel);
        } else if (command.equals(showEnergy)) {
            skyPanel.displayEnergy(true);
            energyButton.setText(hideEnergy);
        } else if (command.equals(showPos)) {
            skyPanel.displayPos(true);
            positionButton.setText(hidePos);
        } else if (command.equals(freeze)) {
            skyPanel.setFreeze(true);
            freezeButton.setText(unFreeze);
        } else if (command.equals(hideAcc)) {
            skyPanel.displayAcc(false);
            accButton.setText(showAcc);
        } else if (command.equals(hideVel)) {
            skyPanel.displayVel(false);
            velButton.setText(showVel);
        } else if (command.equals(hideEnergy)) {
            skyPanel.displayEnergy(false);
            energyButton.setText(showEnergy);
        } else if (command.equals(hidePos)) {
            skyPanel.displayPos(false);
            skyPanel.resetPos();
            positionButton.setText(showPos);
        } else if (command.equals(unFreeze)) {
            skyPanel.setFreeze(false);
            freezeButton.setText(freeze);
        }
    }


    public SkyPanel getPanel() {
        return skyPanel;
    }


    public void stateChanged(ChangeEvent event) {
        JSlider source = (JSlider) event.getSource();
        if (!source.getValueIsAdjusting()) sun.setMass(slider.getValue());
    }
}


class SkyPanel extends JPanel {

    private final int numberOfPlanets;

    private final int maxNumberOfOrbits;

    private final int xMaxPixels;

    private final double maxPlanetDia = 60;

    private final boolean[] orbitOccupied;

    private final Planet[] planet;

    private final Vector<Stars> starGroup;

    private final Sun sun;

    private boolean displayEnergy = false, displayVel = false, displayAcc = false, displayPos = false, freeze = false;


    public SkyPanel(double dt, int mMaxXPosition, int mNumberOfPlanets, Sun mSun) {
        xMaxPixels = mMaxXPosition;
        numberOfPlanets = mNumberOfPlanets;
        sun = mSun;

        int numVisibleOrbits = (int) (((double) xMaxPixels / sun.getAbsPos(Coord.X) - sun.getDiameter() / 2) / maxPlanetDia);

        maxNumberOfOrbits = Math.max(numVisibleOrbits, numberOfPlanets);

        orbitOccupied = new boolean[maxNumberOfOrbits];
        for (int i = 0; i < maxNumberOfOrbits; ++i)
            orbitOccupied[i] = false;

        planet = new Planet[numberOfPlanets];

        for (int i = 0; i < numberOfPlanets; ++i) {
            double minPlanetDia = 20;
            double diameter = minPlanetDia + Math.random() * (maxPlanetDia - minPlanetDia);
            planet[i] = new Planet(dt, // Time increment for calculating motion.
                    diameter, // Diameter of planet
                    getUnoccupiedOrbit(), // Radius of orbit for planet
                    SolarSystem.getColor(i), // Color of planet
                    diameter / maxPlanetDia, // Mass of planet
                    this, // Host panel
                    sun);
        }

        starGroup = new Vector<>(1);
        populateStarGroup(1);
    }


    public void populateStarGroup(int mLengthScaleFactor) {

        if (starGroup.size() >= mLengthScaleFactor) return;


        int numStars = 30;
        Stars stars = new Stars(numStars / SkyFrame.getScaleFactor(), // number of stars
                SkyFrame.getScaleFactor() * (xMaxPixels - sun.getAbsPos(Coord.X)), // maximum distance from sun
                (SkyFrame.getScaleFactor() - 1) * (xMaxPixels - sun.getAbsPos(Coord.X)), // minimum distance from sun

                sun);
        starGroup.addElement(stars);
    }


    private double getUnoccupiedOrbit() {
        int trialOrbitNum = (int) (Math.random() * maxNumberOfOrbits);

        while (true) {
            if (!orbitOccupied[trialOrbitNum]) {
                orbitOccupied[trialOrbitNum] = true;
                break;
            } else trialOrbitNum = (int) (Math.random() * maxNumberOfOrbits);
        }

        return 1.5 * sun.getDiameter() + trialOrbitNum * maxPlanetDia;
    }


    public void paintComponent(Graphics comp) {
        super.paintComponent(comp);
        Graphics2D comp2D = (Graphics2D) comp;

        sun.draw(comp2D);

        for (int i = 0; i < SkyFrame.getScaleFactor(); ++i)
            (starGroup.get(i)).draw(comp2D);

        for (int i = 0; i < numberOfPlanets; ++i) {
            if (!freeze) planet[i].translate();
            planet[i].draw(comp2D);
            if (displayVel) drawVelocityVector(comp2D, i);
            if (displayAcc) drawAccelerationVector(comp2D, i);
            if (displayEnergy) drawEnergy(comp2D, i);
        }
    }


    private void drawEnergy(Graphics2D comp2D, int i) {
        comp2D.setColor(Color.white);


        String percentKE = Integer.toString((int) (100 * planet[i].getKE() / (planet[i].getKE() + planet[i].getPE())));
        comp2D.drawString("KE% " + percentKE, (int) planet[i].getAbsPos(Coord.X) + 50, (int) planet[i].getAbsPos(Coord.Y) + 70);

        String percentPE = Integer.toString((int) (100 * planet[i].getPE() / (planet[i].getKE() + planet[i].getPE())));
        comp2D.drawString("PE% " + percentPE, (int) planet[i].getAbsPos(Coord.X) + 50, (int) planet[i].getAbsPos(Coord.Y) + 82);
    }


    private void drawVelocityVector(Graphics2D comp2D, int i) {
        float xPos = (float) planet[i].getAbsPos(Coord.X);
        float yPos = (float) planet[i].getAbsPos(Coord.Y);
        float xVel = xPos + 25 * (float) planet[i].getVel(Coord.X) / (float) SkyFrame.getScaleFactor();
        float yVel = yPos + 25 * (float) planet[i].getVel(Coord.Y) / (float) SkyFrame.getScaleFactor();

        comp2D.setColor(Color.white);
        comp2D.draw(new Line2D.Float(xPos, yPos, xVel, yVel));
        drawArrow(comp2D, xPos, yPos, xVel, yVel);
    }

    private void drawAccelerationVector(Graphics2D comp2D, int i) {
        double xPos = planet[i].getAbsPos(Coord.X);
        double yPos = planet[i].getAbsPos(Coord.Y);
        double xAcc = xPos + 2000 * planet[i].getAcc(Coord.X) / SkyFrame.getScaleFactor();
        double yAcc = yPos + 2000 * planet[i].getAcc(Coord.Y) / SkyFrame.getScaleFactor();

        comp2D.setColor(Color.white);
        comp2D.draw(new Line2D.Float((float) xPos, (float) yPos, (float) xAcc, (float) yAcc));
        drawArrow(comp2D, (float) xPos, (float) yPos, (float) xAcc, (float) yAcc);
    }


    private void drawArrow(Graphics2D comp2D, float x0, float y0, float x1, float y1) {
        float dx = x1 - x0, dy = y1 - y0;
        float r = (float) Math.sqrt(dx * dx + dy * dy);
        float len = 15 / SkyFrame.getScaleFactor();
        int[] x, y;
        int numPts = 4;

        x = new int[numPts];
        y = new int[numPts];
        x[0] = (int) x1;
        y[0] = (int) y1;
        x[1] = (int) (x1 - len * dy / (2 * r));
        y[1] = (int) (y1 + len * dx / (2 * r));
        x[2] = (int) (x1 + len * dx / r);
        y[2] = (int) (y1 + len * dy / r);
        x[3] = (int) (x1 + len * dy / (2 * r));
        y[3] = (int) (y1 - len * dx / (2 * r));

        GeneralPath path = new GeneralPath(GeneralPath.WIND_EVEN_ODD, numPts);
        path.moveTo(x[0], y[0]);
        for (int i = 1; i < numPts; ++i)
            path.lineTo(x[i], y[i]);
        path.closePath();
        comp2D.fill(path);
    }


    public void displayAcc(boolean value) {
        displayAcc = value;
    }


    public void displayVel(boolean value) {
        displayVel = value;
    }


    public void displayEnergy(boolean value) {
        displayEnergy = value;
    }


    public void displayPos(boolean value) {
        displayPos = value;
    }


    public boolean showPos() {
        return displayPos;
    }


    public void resetPos() {
        for (int i = 0; i < numberOfPlanets; ++i)
            planet[i].erasePrevPositions();
    }


    public void setFreeze(boolean value) {
        freeze = value;
    }


}


class Stars {

    private final double[] xPos;
    private final double[] yPos;

    private final int numOfStars;

    private final Sun sun;


    public Stars(int mNumOfStars, double maxRadius, double minRadius, Sun mSun) {
        numOfStars = Math.abs(mNumOfStars);
        xPos = new double[numOfStars];
        yPos = new double[numOfStars];
        maxRadius = Math.abs(maxRadius);
        minRadius = Math.abs(minRadius);
        sun = mSun;

        if (minRadius == 0) minRadius += sun.getDiameter();

        for (int i = 0; i < numOfStars; ++i) {
            xPos[i] = (float) Math.random() * maxRadius;
            float yMin = xPos[i] > minRadius ? 0 : (float) Math.sqrt(minRadius * minRadius - xPos[i] * xPos[i]);
            float yMax = (float) Math.sqrt(maxRadius * maxRadius - xPos[i] * xPos[i]);
            yPos[i] = (float) Math.random() * (yMax - yMin) + yMin;
            if (Math.random() > 0.5) xPos[i] *= -1;
            if (Math.random() > 0.5) yPos[i] *= -1;
        }
    }


    public void draw(Graphics2D comp2D) {
        comp2D.setColor(Color.yellow);


        for (int i = 0; i < numOfStars; ++i) {
            double diameter = 5;
            double actualDia = diameter * Math.random();

            Ellipse2D.Double starShape = new Ellipse2D.Double(xPos[i] / SkyFrame.getScaleFactor() + sun.getAbsPos(Coord.X), yPos[i] / SkyFrame.getScaleFactor() + sun.getAbsPos(Coord.Y), actualDia, actualDia);
            comp2D.fill(starShape);
        }
    }
}


abstract class Satellite implements Host {

    protected final Color color;

    protected double xPos, yPos;

    protected double xVel = 0, yVel = 0;

    protected double xAcc = 0, yAcc = 0, dt;

    protected double diameter;

    protected double mass;

    protected Host host;


    public Satellite(double mDt, double mDiameter, double mRadius, Color mColor, double mMass, Host mHost) {
        dt = mDt;
        diameter = mDiameter;
        mass = mMass;
        color = mColor;
        host = mHost;

        xPos = Math.random() * mRadius;
        yPos = Math.sqrt(mRadius * mRadius - xPos * xPos);

        if (probOneHalf()) xPos *= -1;
        if (probOneHalf()) yPos *= -1;
    }


    protected boolean probOneHalf() {
        return Math.random() > 0.5;
    }


    public double getAbsPos(byte mAxis) {
        if (mAxis == Coord.X) return xPos / SkyFrame.getScaleFactor() + host.getAbsPos(mAxis);
        else return yPos / SkyFrame.getScaleFactor() + host.getAbsPos(mAxis);
    }


    protected double getRelPosUnscaled(byte mAxis) {
        if (mAxis == Coord.X) return xPos;
        else return yPos;
    }


    public double getTheta() {
        double theta = Math.atan(yPos / xPos);
        if (xPos < 0) theta = Math.PI + theta;
        else if (xPos > 0 && yPos < 0) theta = 2 * Math.PI + theta;
        return theta;
    }

    public double getVel(byte mAxis) {
        if (mAxis == Coord.X) return xVel;
        else return yVel;
    }

    public double getAcc(byte mAxis) {
        if (mAxis == Coord.X) return xAcc;
        else return yAcc;
    }


    public double getDiameter() {
        return diameter;
    }


    public double getMass() {
        return mass;
    }

    public void draw(Graphics2D comp2D) {
        double zoomedDiameter = diameter / SkyFrame.getScaleFactor();
        Ellipse2D.Double form = new Ellipse2D.Double((int) (getAbsPos(Coord.X) - zoomedDiameter / 2), (int) (getAbsPos(Coord.Y) - zoomedDiameter / 2), (int) zoomedDiameter, (int) zoomedDiameter);
        comp2D.setColor(color);
        comp2D.fill(form);
    }
}


class Planet extends Satellite {

    private final int maxPrevPositions = 128;

    private final byte numOfMoons;

    private final int[] prevXPos;
    private final int[] prevYPos;

    private final SkyPanel hostPanel;

    private int currPosition = 0;

    private int count = 0;

    private boolean looped = false;

    private Moon[] moon;


    public Planet(double mDt, double mDiameter, double mRadius, Color mColor, double mMass, SkyPanel mHostPanel, Host mHost) {
        super(mDt, mDiameter, mRadius, mColor, mMass, mHost);
        hostPanel = mHostPanel;
        prevXPos = new int[maxPrevPositions];
        prevYPos = new int[maxPrevPositions];

        for (int i = 0; i < maxPrevPositions; ++i)
            prevXPos[i] = prevYPos[i] = 0;

        byte maxNumOfMoons = 4;
        numOfMoons = (byte) (Math.random() * maxNumOfMoons);

        if (numOfMoons > 0) {
            moon = new Moon[numOfMoons];
            for (int i = 0; i < numOfMoons; ++i)
                moon[i] = new Moon(dt, // Time increment for calculating motion.
                        Math.random() * 6 + 3, // Moon's diameter.
                        diameter * (1 + 0.2 * Math.random()),// Radius of orbit.
                        this); // Host about which moon orbits.
        }

        double distToHost = Math.sqrt(getRelPosUnscaled(Coord.X) * getRelPosUnscaled(Coord.X) + getRelPosUnscaled(Coord.Y) * getRelPosUnscaled(Coord.Y));

        double forceDueToGravity = host.getMass() * mass / (distToHost * distToHost);
        double theta = getTheta();

        yVel = Math.sqrt(host.getMass() / distToHost) * Math.sin(theta + Math.PI / 2);
        xVel = Math.sqrt(host.getMass() / distToHost) * Math.cos(theta + Math.PI / 2);

        xAcc = forceDueToGravity * Math.cos(theta + Math.PI) / mass;
        yAcc = forceDueToGravity * Math.sin(theta + Math.PI) / mass;

        double TE = mass * (xVel * xVel + yVel * yVel) / 2 + mass * Math.sqrt(xAcc * xAcc + yAcc * yAcc) * distToHost;
    }


    public void translate() {
        double force, distToHost, theta;

        for (int i = 0; i < 1 / dt; ++i) {
            distToHost = Math.sqrt(getRelPosUnscaled(Coord.X) * getRelPosUnscaled(Coord.X) + getRelPosUnscaled(Coord.Y) * getRelPosUnscaled(Coord.Y));

            force = mass * host.getMass() / (distToHost * distToHost);
            theta = getTheta();

            xAcc = force * Math.cos(theta + Math.PI) / mass;
            yAcc = force * Math.sin(theta + Math.PI) / mass;

            xPos += xVel * dt + dt * dt * xAcc / 2;
            yPos += yVel * dt + dt * dt * yAcc / 2;

            xVel += xAcc * dt;
            yVel += yAcc * dt;

            if (hostPanel.showPos()) recordPosition();

            for (int j = 0; j < numOfMoons; ++j)
                moon[j].translate();
        }
    }


    private void recordPosition() {
        int recordPosition = 5;
        if (++count > (int) (recordPosition / dt)) {
            count = 0;
            if (currPosition == maxPrevPositions) {
                currPosition = 0;
                looped = true;
            }

            prevXPos[currPosition] = (int) getRelPosUnscaled(Coord.X);
            prevYPos[currPosition++] = (int) getRelPosUnscaled(Coord.Y);
        }
    }


    public double getKE() {
        return mass * (xVel * xVel + yVel * yVel) / 2;
    }


    public double getPE() {
        return mass * Math.sqrt(xAcc * xAcc + yAcc * yAcc) * Math.sqrt(xPos * xPos + yPos * yPos);
    }


    public void draw(Graphics2D comp2D) {
        double zoomedDiameter = diameter / SkyFrame.getScaleFactor();

        comp2D.setColor(Color.darkGray);
        comp2D.fill(semiCircle(getAbsPos(Coord.X), getAbsPos(Coord.Y), zoomedDiameter, getTheta() + 3 * Math.PI / 2));

        comp2D.setColor(color);
        comp2D.fill(semiCircle(getAbsPos(Coord.X), getAbsPos(Coord.Y), zoomedDiameter, getTheta() + Math.PI / 2));

        if (hostPanel.showPos()) drawPrevPositions(comp2D);

        for (int i = 0; i < numOfMoons; ++i)
            moon[i].draw(comp2D);
    }


    private void drawPrevPositions(Graphics2D comp2D) {
        int end;
        int x, y;
        if (looped) end = maxPrevPositions;
        else end = currPosition;

        for (int i = 0; i < end; ++i) {
            x = prevXPos[i] / SkyFrame.getScaleFactor() + (int) host.getAbsPos(Coord.X);
            y = prevYPos[i] / SkyFrame.getScaleFactor() + (int) host.getAbsPos(Coord.Y);
            comp2D.fill(new Ellipse2D.Float(x, y, 1, 1));
        }
    }


    public void erasePrevPositions() {
        looped = false;
        currPosition = 0;
        count = 0;
    }


    private Area semiCircle(double x, double y, double dia, double theta) {
        double xp = x - dia / 2;
        double yp = y - dia / 2;
        Ellipse2D.Double circle = new Ellipse2D.Double(xp, yp, dia, dia);
        Rectangle2D.Double rect = new Rectangle2D.Double(xp, yp, dia, dia / 2);
        Area semiCircleArea = new Area(circle);
        Area rectArea = new Area(rect);
        semiCircleArea.subtract(rectArea);
        semiCircleArea.transform(AffineTransform.getRotateInstance(theta, x, y));
        return semiCircleArea;
    }
}


class Moon extends Satellite {

    public Moon(double mDt, double mDiameter, double mRadius, Host mHost) {
        super(mDt, mDiameter, mRadius, Color.lightGray, 1, mHost);
    }


    public void translate() {
        double r, theta;

        r = Math.sqrt(getRelPosUnscaled(Coord.X) * getRelPosUnscaled(Coord.X) + getRelPosUnscaled(Coord.Y) * getRelPosUnscaled(Coord.Y));


        theta = getTheta() + Math.sqrt(1000 * host.getMass() / (r * r * r)) * dt;

        xPos = r * Math.cos(theta);
        yPos = r * Math.sin(theta);
    }
}


class Sun extends Satellite {

    public Sun(double mDiameter, double mMass, double mXPos, double mYPos) {
        super(1, mDiameter, 0, Color.yellow, mMass, null);

        xPos = mXPos;
        yPos = mYPos;
    }

    public double getAbsPos(byte mAxis) {
        if (mAxis == Coord.X) return xPos;
        else return yPos;
    }


    public void setMass(double newMass) {
        diameter *= newMass / mass;
        mass = newMass;
    }

    public void draw(Graphics2D comp2D) {
        super.draw(comp2D);

        BasicStroke pen = new BasicStroke(2F);
        comp2D.setStroke(pen);

        double xf, yf;

        double zoomedDiameter = diameter / SkyFrame.getScaleFactor();

        for (int i = 0; i < 10; ++i) {
            if (probOneHalf()) xf = xPos + Math.random() * zoomedDiameter * 1.5;
            else xf = xPos - Math.random() * zoomedDiameter * 1.5;

            if (probOneHalf()) yf = yPos + Math.random() * zoomedDiameter * 1.5;
            else yf = yPos - Math.random() * zoomedDiameter * 1.5;

            Line2D.Double ray = new Line2D.Double(xPos, yPos, xf, yf);
            comp2D.draw(ray);
        }
    }

}
