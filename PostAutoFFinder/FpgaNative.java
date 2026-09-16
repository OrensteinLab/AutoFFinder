package PostAutoFFinder;

import java.nio.ByteBuffer;

final class FpgaNative {
    private static boolean loaded;

    private FpgaNative() {
    }

    static synchronized void loadLibrary(String libraryPath) {
        if (loaded) {
            return;
        }
        if (libraryPath == null || libraryPath.isEmpty()) {
            System.loadLibrary("relev_jni");
        } else {
            System.load(libraryPath);
        }
        loaded = true;
    }

    static native int run(
            String xclbinPath,
            String chromosomePath,
            String patternsPath,
            int editDistance,
            ByteBuffer output);
}