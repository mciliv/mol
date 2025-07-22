import 'package:flutter/material.dart';
import 'package:provider/provider.dart';
import 'camera_service.dart';
import 'analysis_service.dart';
import 'payment_service.dart';

class CameraView extends StatelessWidget {
  @override
  Widget build(BuildContext context) {
    return Consumer<CameraService>(
      builder: (context, cameraService, child) {
        if (!cameraService.isInitialized) {
          return CameraInitializing();
        }
        
        if (cameraService.hasError) {
          return CameraError(error: cameraService.errorMessage);
        }
        
        return CameraPreviewWidget();
      },
    );
  }
}

class CameraInitializing extends StatelessWidget {
  @override
  Widget build(BuildContext context) {
    return Container(
      decoration: BoxDecoration(
        color: Colors.black,
        borderRadius: BorderRadius.circular(8),
      ),
      child: Center(
        child: Column(
          mainAxisAlignment: MainAxisAlignment.center,
          children: [
            CircularProgressIndicator(
              valueColor: AlwaysStoppedAnimation<Color>(Colors.white.withOpacity(0.7)),
            ),
            SizedBox(height: 16),
            Text(
              'Initializing camera...',
              style: TextStyle(
                color: Colors.white.withOpacity(0.7),
                fontSize: 16,
              ),
            ),
          ],
        ),
      ),
    );
  }
}

class CameraError extends StatelessWidget {
  final String error;
  
  const CameraError({required this.error});
  
  @override
  Widget build(BuildContext context) {
    return Container(
      decoration: BoxDecoration(
        color: Colors.black,
        borderRadius: BorderRadius.circular(8),
      ),
      child: Center(
        child: Column(
          mainAxisAlignment: MainAxisAlignment.center,
          children: [
            Icon(
              Icons.camera_alt_outlined,
              size: 64,
              color: Colors.white.withOpacity(0.3),
            ),
            SizedBox(height: 16),
            Text(
              'Camera Error',
              style: TextStyle(
                color: Colors.white,
                fontSize: 18,
                fontWeight: FontWeight.w500,
              ),
            ),
            SizedBox(height: 8),
            Text(
              error,
              style: TextStyle(
                color: Colors.white.withOpacity(0.7),
                fontSize: 14,
              ),
              textAlign: TextAlign.center,
            ),
            SizedBox(height: 16),
            ElevatedButton(
              onPressed: () {
                context.read<CameraService>().initializeCamera();
              },
              child: Text('Retry'),
            ),
          ],
        ),
      ),
    );
  }
}

class CameraPreviewWidget extends StatelessWidget {
  @override
  Widget build(BuildContext context) {
    return Consumer2<CameraService, AnalysisService>(
      builder: (context, cameraService, analysisService, child) {
        return Container(
          decoration: BoxDecoration(
            color: Colors.black,
            borderRadius: BorderRadius.circular(8),
          ),
          child: Stack(
            children: [
              // Camera preview - MUCH SIMPLER THAN WEB!
              ClipRRect(
                borderRadius: BorderRadius.circular(8),
                child: AspectRatio(
                  aspectRatio: 16 / 9,
                  child: Container(
                    color: Colors.grey[800], // Placeholder for camera preview
                    child: Center(
                      child: Text(
                        'Camera Preview\n(CameraPreview widget)',
                        style: TextStyle(
                          color: Colors.white.withOpacity(0.7),
                          fontSize: 16,
                        ),
                        textAlign: TextAlign.center,
                      ),
                    ),
                  ),
                ),
              ),
              
              // Crosshair overlay
              Positioned.fill(
                child: Center(
                  child: Container(
                    width: 60,
                    height: 60,
                    decoration: BoxDecoration(
                      border: Border.all(
                        color: Colors.white.withOpacity(0.7),
                        width: 1, // Minimal border per ui.mdc
                      ),
                      shape: BoxShape.circle,
                    ),
                  ),
                ),
              ),
              
              // Capture instruction
              Positioned(
                bottom: 20,
                left: 0,
                right: 0,
                child: Center(
                  child: Container(
                    padding: EdgeInsets.symmetric(horizontal: 20, vertical: 12),
                    decoration: BoxDecoration(
                      color: Colors.black.withOpacity(0.85), // Overlay per ui.mdc
                      borderRadius: BorderRadius.circular(8),
                    ),
                    child: Text(
                      'Center and tap',
                      style: TextStyle(
                        color: Colors.white,
                        fontSize: 14,
                      ),
                      textAlign: TextAlign.center,
                    ),
                  ),
                ),
              ),
              
              // Capture button (invisible overlay for tap)
              Positioned.fill(
                child: GestureDetector(
                  onTap: () => _captureImage(context, cameraService, analysisService),
                  child: Container(color: Colors.transparent),
                ),
              ),
              
              // Switch camera button (if multiple cameras available)
              if (cameraService.hasMultipleCameras)
                Positioned(
                  top: 16,
                  right: 16,
                  child: IconButton(
                    icon: Icon(Icons.flip_camera_android, color: Colors.white),
                    onPressed: () => cameraService.switchCamera(),
                    style: IconButton.styleFrom(
                      backgroundColor: Colors.black.withOpacity(0.5),
                      minimumSize: Size(40, 40), // Touch-friendly
                    ),
                  ),
                ),
            ],
          ),
        );
      },
    );
  }
  
  Future<void> _captureImage(
    BuildContext context,
    CameraService cameraService,
    AnalysisService analysisService,
  ) async {
    try {
      final imageFile = await cameraService.captureImage();
      if (imageFile != null) {
        await analysisService.analyzeImage(imageFile);
      }
    } catch (error) {
      ScaffoldMessenger.of(context).showSnackBar(
        SnackBar(content: Text('Capture failed: $error')),
      );
    }
  }
} 