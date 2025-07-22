import 'package:flutter/foundation.dart';
import 'package:camera/camera.dart';
import 'dart:io';

class CameraService extends ChangeNotifier {
  CameraController? _controller;
  List<CameraDescription> _cameras = [];
  bool _isInitialized = false;
  bool _hasError = false;
  String _errorMessage = '';
  int _currentCameraIndex = 0;
  
  // Getters
  bool get isInitialized => _isInitialized;
  bool get hasError => _hasError;
  String get errorMessage => _errorMessage;
  bool get hasMultipleCameras => _cameras.length > 1;
  CameraController? get controller => _controller;
  
  CameraService() {
    initializeCamera();
  }
  
  // Initialize camera
  Future<void> initializeCamera() async {
    try {
      _hasError = false;
      _errorMessage = '';
      notifyListeners();
      
      // Get available cameras
      _cameras = await availableCameras();
      
      if (_cameras.isEmpty) {
        throw Exception('No cameras available on this device');
      }
      
      // Use back camera by default, fallback to first available
      _currentCameraIndex = _cameras.indexWhere(
        (camera) => camera.lensDirection == CameraLensDirection.back
      );
      if (_currentCameraIndex == -1) _currentCameraIndex = 0;
      
      // Initialize controller
      _controller = CameraController(
        _cameras[_currentCameraIndex],
        ResolutionPreset.high,
        enableAudio: false, // No audio needed for molecular analysis
      );
      
      await _controller!.initialize();
      _isInitialized = true;
      notifyListeners();
      
    } catch (error) {
      _hasError = true;
      _errorMessage = error.toString();
      _isInitialized = false;
      notifyListeners();
    }
  }
  
  // Switch between front and back cameras
  Future<void> switchCamera() async {
    if (_cameras.length <= 1) return;
    
    try {
      await _controller?.dispose();
      _currentCameraIndex = (_currentCameraIndex + 1) % _cameras.length;
      
      _controller = CameraController(
        _cameras[_currentCameraIndex],
        ResolutionPreset.high,
        enableAudio: false,
      );
      
      await _controller!.initialize();
      notifyListeners();
      
    } catch (error) {
      _hasError = true;
      _errorMessage = 'Failed to switch camera: $error';
      notifyListeners();
    }
  }
  
  // Capture image
  Future<File?> captureImage() async {
    if (!_isInitialized || _controller == null) {
      throw Exception('Camera not initialized');
    }
    
    try {
      final XFile image = await _controller!.takePicture();
      return File(image.path);
    } catch (error) {
      throw Exception('Failed to capture image: $error');
    }
  }
  
  @override
  void dispose() {
    _controller?.dispose();
    super.dispose();
  }
} 