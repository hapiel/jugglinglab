// OpenFilesServer.java
//
// Copyright 2002-2021 Jack Boyce and the Juggling Lab contributors

package jugglinglab.util;

import java.io.*;
import java.net.*;
import java.text.MessageFormat;
import java.util.ResourceBundle;

import jugglinglab.core.ApplicationWindow;
import jugglinglab.core.Constants;


public class OpenFilesServer extends Thread {
    static final ResourceBundle guistrings = jugglinglab.JugglingLab.guistrings;
    static final ResourceBundle errorstrings = jugglinglab.JugglingLab.errorstrings;

    static protected final int OPEN_FILES_PORT = 8686;

    protected ServerSocket listen_socket;


    public OpenFilesServer() {
        try {
            listen_socket = new ServerSocket(OPEN_FILES_PORT);
        } catch (IOException e) {
            if (Constants.DEBUG_OPEN_SERVER)
                System.out.println("Server already running on machine; thread is not starting");
            // ErrorDialog.handleFatalException(e);
            return;
        }
        start();
    }

    // Server thread loops forever, listening for connections on our port
    public void run() {
        if (Constants.DEBUG_OPEN_SERVER)
            System.out.println("Server: listening on port " + OPEN_FILES_PORT);

        try {
            while (true) {
                Socket client_socket = listen_socket.accept();

                if (Constants.DEBUG_OPEN_SERVER) {
                    System.out.println("Server got a connection from " +
                            client_socket.getInetAddress() + ":" +
                            client_socket.getPort());
                }

                if (client_socket.getInetAddress().toString().contains("127.0.0.1"))
                    new Connection(client_socket);
                else if (Constants.DEBUG_OPEN_SERVER) {
                    System.out.println("Ignoring connection request; not from 127.0.0.1");
                }
            }
        } catch (IOException e) {
            ErrorDialog.handleFatalException(e);
        }
    }

    // Try to signal another instance of Juggling Lab on this machine to open
    // the file. If the open command is successfully handed off, return true.
    // Otherwise return false.
    static public boolean tryOpenFile(File f) {
        Socket s = null;
        BufferedReader sin = null;
        PrintStream sout = null;

        try {
            s = new Socket("localhost", OPEN_FILES_PORT);
            sin = new BufferedReader(new InputStreamReader(s.getInputStream()));
            sout = new PrintStream(s.getOutputStream());

            if (Constants.DEBUG_OPEN_SERVER) {
                System.out.println("Connected to " + s.getInetAddress()
                           + ":"+ s.getPort());
            }

            String line;

            sout.println("identify");
            line = sin.readLine();
            if (!line.equals("Juggling Lab version " + Constants.version)) {
                if (Constants.DEBUG_OPEN_SERVER) {
                    System.out.println("ID response didn't match: " + line);
                    System.out.println("exiting");
                }
                return false;
            }

            sout.println("open " + f);
            line = sin.readLine();
            if (!line.startsWith("opening ")) {
                if (Constants.DEBUG_OPEN_SERVER) {
                    System.out.println("Open response didn't match: " + line);
                    System.out.println("exiting");
                }
                return false;
            }

            sout.println("done");
            line = sin.readLine();

            return true;

            /*
            // interactive console app to talk with the server

            BufferedReader in = new BufferedReader(new InputStreamReader(System.in));
            while (true) {
                System.out.print("> ");
                System.out.flush();

                line = in.readLine();
                if (line == null)
                    break;

                sout.println(line);
                line = sin.readLine();
                if (line == null) {
                    System.out.println("Connection closed by server.");
                    break;
                }
                System.out.println(line);
            }
            return true;
            */
        } catch (IOException e) {
            if (Constants.DEBUG_OPEN_SERVER)
                System.out.println("No server detected on port " + OPEN_FILES_PORT);
            return false;
        } finally {
            try {
                if (s != null) s.close();
            } catch (IOException e2) {}
        }
    }
}

// Server thread that handles all communication with a client.

class Connection extends Thread {
    static final ResourceBundle guistrings = jugglinglab.JugglingLab.guistrings;
    static final ResourceBundle errorstrings = jugglinglab.JugglingLab.errorstrings;

    protected Socket client;
    protected BufferedReader in;
    protected PrintStream out;


    public Connection(Socket client_socket) {
        client = client_socket;
        try {
            in = new BufferedReader(new InputStreamReader(client.getInputStream()));
            out = new PrintStream(client.getOutputStream());
        } catch (IOException ioe) {
            try {
                client.close();
            } catch (IOException ioe2) {}
            ErrorDialog.handleFatalException(ioe);
            return;
        }
        start();
    }

    public void run() {
        String line;
        StringBuffer revline;
        int len;

        if (Constants.DEBUG_OPEN_SERVER)
            System.out.println("Server started connection thread");

        try {
            while (true) {
                line = in.readLine();
                if (line == null)
                    return;

                if (line.startsWith("open ")) {
                    String filepath = line.substring(5, line.length());
                    File file = new File(filepath);

                    out.println("opening " + filepath);

                    try {
                        ApplicationWindow.openJMLFile(file);
                    } catch (JuggleExceptionUser jeu) {
                        String template = errorstrings.getString("Error_reading_file");
                        Object[] arguments = { file.getName() };
                        String msg = MessageFormat.format(template, arguments) +
                                     ":\n" + jeu.getMessage();
                        new ErrorDialog(null, msg);
                    } catch (JuggleExceptionInternal jei) {
                        ErrorDialog.handleFatalException(jei);
                    }
                } else if (line.startsWith("identify")) {
                    out.println("Juggling Lab version " + Constants.version);
                } else if (line.startsWith("done")) {
                    out.println("goodbye");
                    return;
                } else
                    out.println(line);
            }
        } catch (IOException ioe) {
        } finally {
            try {
                client.close();
            } catch (IOException ioe2) {}
        }
    }
}